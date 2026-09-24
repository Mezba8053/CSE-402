"""Command-line entry point for the vortex-flowmeter project."""

from pathlib import Path
<<<<<<< HEAD
=======
import os
>>>>>>> 56529ac (added pylbm model for optimized side radius)
import sys


local_dependencies = Path(__file__).resolve().parent / ".python_packages"
if local_dependencies.is_dir():
    sys.path.insert(0, str(local_dependencies))
<<<<<<< HEAD
=======
    existing_pythonpath = os.environ.get("PYTHONPATH", "")
    os.environ["PYTHONPATH"] = str(local_dependencies) + (
        os.pathsep + existing_pythonpath if existing_pythonpath else ""
    )
>>>>>>> 56529ac (added pylbm model for optimized side radius)

from vortexflow.cli import main


if __name__ == "__main__":
    main()
