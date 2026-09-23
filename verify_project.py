"""Run the complete automated verification suite with local dependencies."""

from pathlib import Path
import sys
import unittest


ROOT = Path(__file__).resolve().parent
LOCAL_DEPENDENCIES = ROOT / ".python_packages"
if LOCAL_DEPENDENCIES.is_dir():
    sys.path.insert(0, str(LOCAL_DEPENDENCIES))
sys.path.insert(0, str(ROOT))


def main() -> int:
    suite = unittest.defaultTestLoader.discover(str(ROOT / "tests"))
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    return 0 if result.wasSuccessful() else 1


if __name__ == "__main__":
    raise SystemExit(main())
