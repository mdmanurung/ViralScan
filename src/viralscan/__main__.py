"""Entry point for ``viralscan`` and ``python -m viralscan``.

Keeps ``--version`` and the Python guard free of the heavy imports in ``menu``.
"""

import sys


def main() -> None:
    if sys.version_info < (3, 10):  # noqa: UP036 (old interpreters ignore requires-python)
        sys.exit(f"viralscan needs Python >= 3.10 (found {sys.version.split()[0]}).")
    if sys.argv[1:2] == ["--version"]:
        from viralscan import __version__

        print(f"viralscan {__version__}")
        return
    from viralscan.menu import main as menu_main

    menu_main()


if __name__ == "__main__":
    main()
