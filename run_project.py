"""Command-line entry point for the vortex-flowmeter project."""

from pathlib import Path
import sys


local_dependencies = Path(__file__).resolve().parent / ".python_packages"
if local_dependencies.is_dir():
    sys.path.insert(0, str(local_dependencies))

from vortexflow.cli import main


if __name__ == "__main__":
    main()
