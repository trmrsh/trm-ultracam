"""
Tests for trm.ultracam
"""

import unittest
import numpy as np
from trm import ultracam

class TestWindow(unittest.TestCase):

    def setUp(self):
        """
        Create a couple of compatible Windows
        """
        llx, lly, xbin, ybin = 12, 22, 2, 3
        dat1  = np.zeros((10,20))
        dat1 += 2.
        self.win1 = ultracam.Window(dat1, llx, lly, xbin, ybin) 
        dat2  = np.zeros((10,20))
        dat2 += 3.
        self.win2 = ultracam.Window(dat2, llx, lly, xbin, ybin) 

    def test_index(self):
        self.assertEqual(self.win1[4,5],2.)

    def test_nx(self):
        self.assertEqual(self.win1.nx, 20)

    def test_ny(self):
        self.assertEqual(self.win1.ny, 10)

    def test_llx(self):
        self.assertEqual(self.win1.llx, 12)

    def test_lly(self):
        self.assertEqual(self.win1.lly, 22)

    def test_xbin(self):
        self.assertEqual(self.win1.xbin, 2)

    def test_xbin(self):
        self.assertEqual(self.win1.ybin, 3)

    def test_size(self):
        self.assertEqual(self.win1.size,200)

    def test_iadd(self):
        self.win1 += 2
        self.assertEqual(self.win1[4,5],4.)

    def test_isub(self):
        self.win1 -= 2
        self.assertEqual(self.win1[4,5],0.)

    def test_imul(self):
        self.win1 *= 2
        self.assertEqual(self.win1[4,5],4.)

    def test_idiv(self):
        self.win1 /= 2.
        self.assertEqual(self.win1[4,5],1.)

    def test_add(self):
        win = self.win1 + self.win2
        self.assertEqual(win[4,5],5.)

    def test_sub(self):
        win = self.win1 - self.win2
        self.assertEqual(win[4,5],-1.)

    def test_mul(self):
        win = self.win1 * self.win2
        self.assertEqual(win[4,5],6.)

    def test_div(self):
        win = self.win1 / self.win2
        self.assertEqual(win[4,5],2./3.)


if __name__ == '__main__':
    unittest.main()
