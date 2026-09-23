"""Expose project-local Python dependencies without modifying global Python."""

from pathlib import Path
import sys


dependency_directory = Path(__file__).resolve().parent / ".python_packages"
if dependency_directory.is_dir():
    sys.path.insert(0, str(dependency_directory))
