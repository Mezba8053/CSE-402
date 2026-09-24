"""Small compatibility helpers for third-party PyLBM.

PyLBM 0.11 constructs boundary argument dictionaries with ``locals()`` after
dynamic ``exec`` assignments.  Python 3.14 no longer exposes those assignments
through ``locals()``.  This shim restores only the missing array references; it
does not change any LBM equation or boundary algorithm.
"""

from __future__ import annotations

import sys


def apply_python314_boundary_fix(pylbm_module) -> None:
    if sys.version_info < (3, 14):
        return
    if getattr(pylbm_module.boundary, "_vortexflow_python314_fix", False):
        return

    def patch(boundary_class) -> None:
        original = boundary_class._get_args

        def compatible_get_args(self, distributions):
            arguments = original(self, distributions)
            arguments.update(
                {f"iload{index}": value for index, value in enumerate(self.iload)}
            )
            return arguments

        boundary_class._get_args = compatible_get_args

    patch(pylbm_module.boundary.BoundaryMethod)
    patch(pylbm_module.boundary.BouzidiBounceBack)
    pylbm_module.boundary._vortexflow_python314_fix = True
