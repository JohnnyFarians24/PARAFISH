"""Controls that no callback acts on are not shown to the user.

They stay in the layout because the Run Pipeline callback reads some of them as State,
but a visible control that changes nothing tells the user something the platform does
not do. Every id below must sit inside a component with ``display: none``; the live
controls next to them must not.
"""

from __future__ import annotations

import unittest

from ui.pages import page_results, page_setup

HIDDEN_IN_SETUP = (
    "setup-vsearch-sintax",  # the cutoff that takes effect is setup-sintax-cutoff
    "setup-species-threshold",
    "setup-workers",
    "setup-mafft-mode",
    "setup-adjust-direction",
)
HIDDEN_IN_RESULTS = ("results-alpha-table", "results-taxonomy-table", "results-tracking-table")
SHOWN_IN_SETUP = ("setup-sintax-cutoff", "setup-vsearch-minsize", "setup-raxml-bootstraps", "setup-threads")


def _path_to(node, target, trail=()):
    """The chain of components from ``node`` down to the one whose id is ``target``."""
    trail = trail + (node,)
    if getattr(node, "id", None) == target:
        return trail
    children = getattr(node, "children", None)
    if children is None:
        return None
    if not isinstance(children, (list, tuple)):
        children = [children]
    for child in children:
        found = _path_to(child, target, trail)
        if found:
            return found
    return None


def _hidden(tree, control_id):
    path = _path_to(tree, control_id)
    if path is None:
        raise AssertionError(f"{control_id} is not in the layout")
    return any((getattr(n, "style", None) or {}).get("display") == "none" for n in path)


class DeadControlsHiddenTests(unittest.TestCase):
    def test_setup_controls_without_effect_are_hidden(self):
        tree = page_setup.layout()
        for control_id in HIDDEN_IN_SETUP:
            with self.subTest(control_id=control_id):
                self.assertTrue(_hidden(tree, control_id))

    def test_live_setup_controls_stay_visible(self):
        tree = page_setup.layout()
        for control_id in SHOWN_IN_SETUP:
            with self.subTest(control_id=control_id):
                self.assertFalse(_hidden(tree, control_id))

    def test_empty_results_sections_are_hidden(self):
        tree = page_results.layout()
        for control_id in HIDDEN_IN_RESULTS:
            with self.subTest(control_id=control_id):
                self.assertTrue(_hidden(tree, control_id))


if __name__ == "__main__":
    unittest.main()
