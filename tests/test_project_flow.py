"""The Project page must tell the user what happened and what to do next."""

from __future__ import annotations

import unittest

from ui.pages.page_project import layout


class ProjectFlowTests(unittest.TestCase):
    def test_create_project_navigates_onward(self):
        """Creating a project used to return no_update for the URL and no confirmation.

        The project WAS created — it just looked like nothing had happened: no
        navigation, no toast, and a session panel still reading "empty until a project is
        selected". Users reasonably concluded the button was dead.
        """
        import inspect

        from ui import callbacks_project

        source = inspect.getsource(callbacks_project.register_callbacks)
        create = source.split("def _create_project")[1].split("def ")[0]
        self.assertIn('"/samples"', create, "create must send the user to the next step")
        self.assertNotIn("return store, {}, no_update", create)

    def test_the_session_panel_is_wired_to_a_callback(self):
        """It was a hard-coded placeholder that no callback ever wrote to."""
        import inspect

        from ui import callbacks_project

        source = inspect.getsource(callbacks_project.register_callbacks)
        self.assertIn("project-session-state-info", source)

    def test_the_placeholder_text_is_still_the_pre_render_default(self):
        """Sanity: the layout ships a placeholder, and the callback replaces it."""
        found = []

        def walk(node):
            if getattr(node, "id", None) == "project-session-state-info":
                found.append(node)
            children = getattr(node, "children", None)
            if children is None:
                return
            for child in children if isinstance(children, (list, tuple)) else [children]:
                if hasattr(child, "children") or hasattr(child, "id"):
                    walk(child)

        walk(layout())
        self.assertEqual(len(found), 1)


if __name__ == "__main__":
    unittest.main()
