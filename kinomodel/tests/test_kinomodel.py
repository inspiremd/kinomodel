"""
Unit and regression test for the kinomodel package.
"""

# Import package, test suite, and other packages as needed
import kinomodel as ki
import unittest
import sys
from unittest.mock import MagicMock

# clean traceback msgs
sys.tracebacklimit = 0


class KinomodelTestCase(unittest.TestCase):

    #def test_kinomodel_installed():
    # make sure kinomodel has been correctly installed and importable
    #self.assertTrue("kinomodel" in sys.modules)

    def test_comma(self):
        # make sure the input is expected
        ki.input = MagicMock(return_value=7)  # non-string input
        with self.assertRaises(ValueError):
            ki.main()
        ki.input = MagicMock(return_value='apple')  # irrelevant string input
        with self.assertRaises(ValueError):
            ki.main()


if __name__ == '__main__':
    unittest.main()
