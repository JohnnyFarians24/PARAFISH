"""The user-facing surfaces are English, and the Downloads table's filter still discriminates.

Two guards from the PT→EN translation (2026-08-03), each pinning a different risk.

`UserFacingLanguageTests` is the anti-regression net: every string literal the pipeline or
the UI can put in front of a user — labels, tooltips, empty states, log messages, error
messages, CSV cell values — must be English. It scans string CONSTANTS via the AST rather
than grepping lines, so a Portuguese word inside a developer comment or a docstring (which
no user ever sees) does not fail the suite, and a Portuguese label cannot hide behind one.

`DownloadsTableWiringTests` guards the one place the translation touched an identifier
rather than a label. The Downloads table's column ids double as its visible headers (they
are `.capitalize()`d into place), and a DataTable `filter_query` compares one of them
against a literal cell value. Renaming the column without the query — or the query without
the cell value — leaves a table that renders perfectly and greys nothing out: looks-alive,
does-nothing, the defect shape this project keeps meeting. This asserts the three ends
still agree by exercising the real callback against a real run.
"""

from __future__ import annotations

import ast
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Modules whose strings can reach a user. `scripts/` is developer tooling and `tests/`
# describe defects in the language they were found in, so both are out of scope.
SCANNED = ("ui", "core", "analysis", "databases", "preprocess", "phylo", "qc")

# Unambiguous Portuguese words: none of these is also an English word, so a hit is a
# finding rather than a false positive. Deliberately not exhaustive — it only has to be
# dense enough that a translated sentence cannot slip back in whole.
#
# TWO WORDS ARE DELIBERATELY ABSENT, both measured as false positives when this test was
# first run: "com" (matches every github.com/*.com URL) and "define" (an English
# imperative — "Define TRIMMOMATIC_JAR_PATH in the environment" is correct English). A
# detector that cries wolf gets its assertion weakened by the next person; these two
# stay out so every hit it does report is real.
PT_WORDS = (
    "não|nao|nenhum|nenhuma|sem|para|dos|das|uma|este|esta|isso|aqui|onde|quando|porque|"
    "corre|correu|correr|falhou|falha|falta|ficheiro|ficheiros|amostra|amostras|dados|fonte|"
    "fontes|estado|nível|nivel|tamanho|número|numero|conjunto|etapa|passo|próximo|proximo|"
    "anterior|projeto|resultado|resultados|leitura|leituras|árvore|arvore|sequência|sequencia|"
    "sequências|sequencias|taxonomia|análise|analise|método|metodo|memória|memoria|reino|filo|"
    "família|familia|género|genero|espécie|especie|disponível|disponivel|"
    "indisponível|indisponivel|apenas|todos|todas|entre|sobre|antes|depois|ainda|"
    "agora|escolhe|verifica|instala|ativa|gerou|produziu|precisa|carregar|guardar"
)
PT_RE = re.compile(rf"\b({PT_WORDS})\b", re.IGNORECASE)


def _literals(path: Path):
    """(lineno, value) for every string constant that is not a docstring."""
    tree = ast.parse(path.read_text(encoding="utf-8"), filename=str(path))
    docstrings = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.ClassDef, ast.FunctionDef, ast.AsyncFunctionDef)):
            body = getattr(node, "body", None)
            if body and isinstance(body[0], ast.Expr) and isinstance(body[0].value, ast.Constant) and isinstance(body[0].value.value, str):
                docstrings.add(id(body[0].value))
    for node in ast.walk(tree):
        if isinstance(node, ast.Constant) and isinstance(node.value, str) and id(node) not in docstrings:
            yield node.lineno, node.value


class UserFacingLanguageTests(unittest.TestCase):
    def test_no_portuguese_reaches_a_user_facing_string(self):
        offenders = []
        for package in SCANNED:
            for path in sorted((ROOT / package).rglob("*.py")):
                if "__pycache__" in path.parts:
                    continue
                for lineno, value in _literals(path):
                    match = PT_RE.search(value)
                    if match:
                        rel = path.relative_to(ROOT)
                        offenders.append(f"{rel}:{lineno}  «{match.group(0)}»  {value[:90]!r}")
        self.assertEqual(
            offenders,
            [],
            "Portuguese found in user-facing strings:\n  " + "\n  ".join(offenders),
        )


class DownloadsTableWiringTests(unittest.TestCase):
    """The Downloads table's column id, its cell value and its filter_query must agree.

    Exercised against `runs/_fixture/current`, which has some artifacts on disk and some
    absent — so the table legitimately contains BOTH an "available" row and a row carrying
    a reason. That is what makes the assertion meaningful: a filter that discriminates
    nothing would pass on a table where every row is the same.
    """

    @classmethod
    def setUpClass(cls):
        from app import create_app

        cls.app = create_app()

    def _exports_component(self):
        run_dir = ROOT / "runs" / "_fixture" / "current"
        if not run_dir.exists():
            self.skipTest("fixture run not present")
        from core.run_registry import load_run_state

        state = load_run_state(run_dir) or {}
        run_config = dict(state.get("run_config") or {})
        run_config.setdefault("run_id", state.get("run_id"))
        run_config["output_dir"] = str(run_dir)

        for key, callback in self.app.callback_map.items():
            if "results-exports" in key:
                # `__wrapped__` is the plain function; the registered one expects Dash's
                # request context, which no unit test has.
                return callback["callback"].__wrapped__(run_config, None, None)
        self.fail("the Downloads callback (results-exports) was not registered")

    @staticmethod
    def _find_table(component):
        """The DataTable inside the returned component tree."""
        stack = [component]
        while stack:
            node = stack.pop()
            if type(node).__name__ == "DataTable":
                return node
            children = getattr(node, "children", None)
            if isinstance(children, (list, tuple)):
                stack.extend(children)
            elif children is not None:
                stack.append(children)
        return None

    def test_the_filter_query_still_matches_the_column_it_greys_out(self):
        table = self._find_table(self._exports_component())
        self.assertIsNotNone(table, "no DataTable rendered in the Downloads section")

        queries = [
            rule["if"]["filter_query"]
            for rule in table.style_data_conditional
            if isinstance(rule.get("if"), dict) and "filter_query" in rule["if"]
        ]
        target = next((q for q in queries if q.startswith("{") and "!=" in q), None)
        self.assertIsNotNone(target, f"the greying-out rule is gone; queries were {queries}")

        column, value = re.match(r"\{(\w+)\}\s*!=\s*\"([^\"]+)\"", target).groups()

        # (1) the query names a column the table actually has
        column_ids = [c["id"] for c in table.columns]
        self.assertIn(column, column_ids, f"filter_query reads {{{column}}}, columns are {column_ids}")

        # (2) every row carries that column, and (3) the compared value is one the rows
        #     really use — otherwise the rule matches everything or nothing in silence.
        values = {row[column] for row in table.data}
        self.assertTrue(all(column in row for row in table.data))
        self.assertIn(value, values, f"no row is ever {value!r}; the rule can only grey out everything")
        self.assertTrue(
            values - {value},
            f"every row is {value!r}; the rule can never fire and the wiring is unproven",
        )
