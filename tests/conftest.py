"""pytest configuration / session-level fixtures.

``menu.py`` imports ``pyfiglet`` at the top level to render the ASCII banner.
When ``pyfiglet`` is not installed in the test environment we stub it out here
so that the rest of the module can still be imported and tested.  This must
happen before any test module imports ``viralscan.menu``.
"""

import sys
from unittest.mock import MagicMock

# Stub pyfiglet if it is not available in the current environment.
try:
    import pyfiglet  # noqa: F401
except ModuleNotFoundError:
    _mock = MagicMock()
    _mock.figlet_format.return_value = "ViralScan"
    sys.modules["pyfiglet"] = _mock
