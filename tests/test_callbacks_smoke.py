"""Every render callback must survive being invoked with empty inputs.

Importing a module does NOT execute its callback bodies, so a `NameError` (a helper used
but never imported) or a `TypeError` (a function called with the wrong signature) hides
until the callback actually fires in the browser — which is exactly how `details` and the
`composition_figure` name-collision reached a user instead of a test. This exercises each
callback the way dash-renderer does on first paint: real function, empty stores.
"""

from __future__ import annotations

import unittest

from app import create_app


class CallbackSmokeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.app = create_app()

    def test_no_output_is_claimed_twice_without_allow_duplicate(self):
        """Two callbacks on one output is a hard Dash error that breaks the whole app —
        it is what took the Project page down when Results re-claimed results-exports."""
        import collections

        seen = collections.Counter()
        for cb in self.app.callback_map.values():
            for output in cb["output"] if isinstance(cb["output"], list) else [cb["output"]]:
                if getattr(output, "allow_duplicate", False):
                    continue
                seen[(output.component_id, output.component_property)] += 1
        dupes = {key: count for key, count in seen.items() if count > 1}
        self.assertEqual(dupes, {}, f"outputs claimed more than once: {dupes}")

    def test_every_server_callback_runs_with_empty_inputs(self):
        """A NameError or TypeError in a callback body must fail HERE, not in a browser."""
        client = self.app.server.test_client()
        deps = client.get("/_dash-dependencies").get_json()

        failures = []
        for cb in deps:
            if cb.get("clientside_function"):
                continue
            inputs = [{"id": i["id"], "property": i["property"], "value": None} for i in cb["inputs"]]
            if not inputs:
                continue
            state = [{"id": s["id"], "property": s["property"], "value": None} for s in cb["state"]]
            outputs = cb["output"]
            multi = outputs.startswith("..")
            body = {
                "output": outputs,
                "outputs": (
                    [{"id": o.split(".")[0], "property": o.split(".")[1]} for o in outputs.strip(".").split("...")]
                    if multi
                    else {"id": outputs.split(".")[0], "property": outputs.split(".")[1]}
                ),
                "inputs": inputs,
                "state": state,
                "changedPropIds": [f"{inputs[0]['id']}.{inputs[0]['property']}"],
            }
            response = client.post("/_dash-update-component", json=body)
            if response.status_code == 500:
                text = response.get_data(as_text=True)
                # PreventUpdate is a legitimate 500-ish no-op in Dash's protocol; a real
                # code error is a NameError / TypeError / AttributeError in the traceback.
                if any(exc in text for exc in ("NameError", "TypeError", "AttributeError", "KeyError", "ImportError")):
                    first = outputs.split("...")[0].strip(".")
                    snippet = next((ln for ln in text.splitlines() if "Error" in ln), "")
                    failures.append(f"{first}: {snippet.strip()[:120]}")

        self.assertEqual(failures, [], "callbacks that crash on empty inputs:\n" + "\n".join(failures))


if __name__ == "__main__":
    unittest.main()
