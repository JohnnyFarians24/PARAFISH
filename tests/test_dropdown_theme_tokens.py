"""The dropdown theme must keep matching the INSTALLED Dash bundle, not just itself.

Two sessions were lost fixing dropdowns against an assumption — that dcc.Dropdown
renders react-select — that was never checked against the real DOM. It is a native
<button class="dash-dropdown"> whose own stylesheet is keyed on --Dash-* design tokens
(defaults are a light theme), so the theme themes it by redefining those tokens.

A CSS-content check ("the selectors exist in the file") cannot catch the failure mode
that actually bit: a Dash upgrade renaming or adding a --Dash-* token, after which the
theme covers a token the bundle no longer uses — or misses one it now does — and the
control silently reverts to light-on-white. So this guard asserts the theme's defined
token set COVERS the set the installed bundle consumes. It reads the installed bundle,
not our own file, so it fails the day install reality drifts from the theme.

Closure of the dropdown work is the user's VISUAL confirmation, not this test — the
lesson of the two lost sessions is that green here never meant the pixels were right.
This only stops a known, invisible drift.
"""

from __future__ import annotations

import os
import re
import unittest
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
THEME = ROOT / "assets" / "parafish_theme.css"


def _dropdown_bundle() -> Path:
    import dash

    return Path(dash.__file__).parent / "dcc" / "async-dropdown.js"


def _tokens_consumed_by_bundle() -> set[str]:
    """Every --Dash-* token the installed dropdown stylesheet references via var(...)."""
    js = _dropdown_bundle().read_text(encoding="utf-8")
    return set(re.findall(r"var\(\s*(--Dash-[A-Za-z-]+)", js))


def _tokens_defined_in_theme() -> set[str]:
    """Every --Dash-* token the theme DEFINES (declaration `--Dash-...:`), not merely
    references."""
    css = THEME.read_text(encoding="utf-8")
    return set(re.findall(r"(--Dash-[A-Za-z-]+)\s*:", css))


class DropdownThemeTokenCoverageTests(unittest.TestCase):
    def test_the_installed_bundle_is_the_native_button_not_react_select(self):
        """Anchors the whole approach: if a future Dash brings react-select back, the
        --Dash-* token strategy is wrong and this states why the guard would need rework."""
        self.assertTrue(_dropdown_bundle().exists(), "async-dropdown.js not found in the installed dash")
        js = _dropdown_bundle().read_text(encoding="utf-8")
        self.assertIn("dash-dropdown", js)
        self.assertNotIn("Select__control", js, "bundle unexpectedly ships react-select — revisit theming")

    def test_theme_covers_every_dash_token_the_bundle_consumes(self):
        consumed = _tokens_consumed_by_bundle()
        defined = _tokens_defined_in_theme()
        self.assertTrue(consumed, "found no --Dash-* tokens in the bundle — parsing broke")
        missing = consumed - defined
        self.assertEqual(
            missing,
            set(),
            "the installed dropdown bundle consumes --Dash-* tokens the theme does not "
            f"define (a Dash upgrade likely renamed/added them): {sorted(missing)}",
        )

    def test_the_accent_token_is_not_left_at_the_bundle_purple(self):
        """--Dash-Fill-Interactive-Strong drives focus/hover; the bundle default is
        #7f4bc4 (purple). It must be mapped to the theme accent, not left purple."""
        css = THEME.read_text(encoding="utf-8")
        m = re.search(r"--Dash-Fill-Interactive-Strong\s*:\s*([^;]+);", css)
        self.assertIsNotNone(m, "theme does not define --Dash-Fill-Interactive-Strong")
        value = m.group(1).strip().lower()
        self.assertNotIn("7f4bc4", value, "accent left at the bundle's purple")
        self.assertIn("--blue", value, "accent should map to the theme's --blue")


if __name__ == "__main__":
    unittest.main()
