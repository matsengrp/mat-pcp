#!/usr/bin/env python3
"""Single definition of the by-path loader for ``workflow/modules/*.py``.

The pipeline scripts run under several conda envs where ``workflow/modules`` is
not on ``sys.path``, so shared modules are imported by absolute path. This is
the one place that logic lives; each script wires it up with a short bootstrap
(see any script's "by-path module loader" block) because the loader is itself
the mechanism for importing shared code and so cannot be imported via itself.
"""

import importlib.util
import os
import sys


def load_workflow_module(module_name, workflow_dir=None):
    """Load ``workflow/modules/<module_name>.py`` by absolute path.

    Fails fast with a clear ``ImportError`` when the module file is missing,
    rather than the cryptic ``AttributeError`` on ``spec.loader`` that a wrong
    path otherwise raises. ``workflow_dir`` defaults to this file's parent's
    parent (the ``workflow/`` directory), so callers may omit it.
    """
    if workflow_dir is None:
        workflow_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    module_path = os.path.join(workflow_dir, "modules", f"{module_name}.py")

    if not os.path.exists(module_path):
        raise ImportError(f"Cannot find module '{module_name}' at {module_path}")

    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module
