#!/usr/bin/env python3
"""Run the ORIGINAL PORTEKrun.py with the in-memory corrections of
portek_runtime_patch.py. Same arguments as PORTEKrun.py.

usage:
    python portek_patched_cli.py find_k PROJECT --strand_specific
    python portek_patched_cli.py find_enriched PROJECT -k 15 -m 2.0

Requirements:
    - portek_runtime_patch.py must sit in the SAME folder as this script.
    - Run with the Python interpreter of the PORT-EK environment
      (e.g. env/bin/python), so that numpy, pandas, biopython, ... are found.

The PORT-EK repository (the folder containing PORTEKrun.py and portek/) is
located in this order:
    1. the PORTEK_DIR environment variable
    2. the current working directory
    3. the parent folder of this script (e.g. PORT-ssEK-main/ for PORT-ssEK-main/script/)
    4. ./PORT-EK-version-2
"""
import os
import sys
import runpy

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))


def is_repo(path):
    """A valid repository contains PORTEKrun.py and the portek/ package."""
    return os.path.isfile(os.path.join(path, "PORTEKrun.py")) and os.path.isdir(
        os.path.join(path, "portek")
    )


def find_repo():
    env_dir = os.environ.get("PORTEK_DIR")
    if env_dir:
        env_dir = os.path.abspath(env_dir)
        if not is_repo(env_dir):
            sys.exit(
                f"PORTEK_DIR is set to '{env_dir}', but PORTEKrun.py and the "
                f"portek/ folder were not found there."
            )
        return env_dir

    candidates = [
        os.getcwd(),
        os.path.dirname(SCRIPT_DIR),
        os.path.abspath("PORT-EK-version-2"),
    ]
    for candidate in candidates:
        if is_repo(candidate):
            return candidate

    searched = "\n  ".join(candidates)
    sys.exit(
        "PORTEKrun.py and the portek/ folder were not found. Searched:\n  "
        f"{searched}\nSet PORTEK_DIR to the PORT-EK-version-2 / PORT-ssEK folder."
    )


def main():
    repo = find_repo()

    # These two lines were previously chained onto the 'if' statement and
    # therefore never executed when the repository WAS found.
    sys.path.insert(0, repo)        # makes 'import portek' work
    sys.path.insert(0, SCRIPT_DIR)  # makes 'import portek_runtime_patch' work

    try:
        import portek  # noqa: F401
    except ImportError as err:
        sys.exit(
            f"Could not import 'portek' from {repo}:\n  {err}\n"
            "If the message names another module (e.g. Bio, numpy, pandas), "
            "run this script with the PORT-EK environment's interpreter, "
            "e.g.  env/bin/python script/portek_patched_cli.py ..."
        )

    try:
        import portek_runtime_patch
    except ImportError as err:
        sys.exit(
            f"Could not import 'portek_runtime_patch':\n  {err}\n"
            f"Place portek_runtime_patch.py in {SCRIPT_DIR} "
            "(the same folder as this script)."
        )

    portek_runtime_patch.apply_findk_patch()

    runner = os.path.join(repo, "PORTEKrun.py")
    sys.argv[0] = runner
    runpy.run_path(runner, run_name="__main__")


if __name__ == "__main__":
    main()
