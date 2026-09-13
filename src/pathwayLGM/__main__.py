"""Support python -m pathwaylgm using the same entry point as the installed CLI."""

from pathwaylgm.main import main

if __name__ == "__main__":
    raise SystemExit(main())
