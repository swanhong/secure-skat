import ctypes
import unittest

import numpy as np

from davies_qfc import Davies


class FakeQFC:
    def __init__(self, cumulative: float, fault: int) -> None:
        self.cumulative = cumulative
        self.fault = fault
        self.called = False

    def __call__(self, *arguments: object) -> None:
        self.called = True
        arguments[-2]._obj.value = self.fault
        arguments[-1]._obj.value = self.cumulative


class DaviesTest(unittest.TestCase):
    def make_davies(self, cumulative: float, fault: int) -> tuple[Davies, FakeQFC]:
        function = FakeQFC(cumulative, fault)
        davies = object.__new__(Davies)
        davies.function = function
        return davies, function

    def test_rank_one_uses_modified_liu_without_qfc(self) -> None:
        davies, function = self.make_davies(-1, 1)

        result = davies.p_value(3.0, np.array([2.0]), 0.125)

        self.assertEqual(result, (0.125, 1))
        self.assertFalse(function.called)

    def test_successful_qfc_returns_davies_p_value(self) -> None:
        davies, _ = self.make_davies(0.75, 0)

        result = davies.p_value(3.0, np.array([2.0, 1.0]), 0.125)

        self.assertEqual(result, (0.25, 1))

    def test_nonzero_ifault_keeps_valid_p_value(self) -> None:
        davies, _ = self.make_davies(0.75, 1)

        result = davies.p_value(3.0, np.array([2.0, 1.0]), 0.125)

        self.assertEqual(result, (0.25, 0))

    def test_invalid_p_value_uses_modified_liu(self) -> None:
        davies, _ = self.make_davies(-1, 1)

        result = davies.p_value(3.0, np.array([2.0, 1.0]), 0.125)

        self.assertEqual(result, (0.125, 0))


if __name__ == "__main__":
    unittest.main()
