"""
Tests for trm.ultracam.

Run with 'python test_ultracam.py'
"""
from __future__ import absolute_import

import unittest
import numpy as np
import ppgplot as pg
from   trm import ultracam

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

    def test_ybin(self):
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

    def test_addc(self):
        win = self.win1 + 10.
        self.assertEqual(win[4,5],12.)

    def test_subc(self):
        win = self.win1 - 10.
        self.assertEqual(win[4,5],-8.)

    def test_mulc(self):
        win = self.win1 * 10.
        self.assertEqual(win[4,5],20.)

    def test_divc(self):
        win = self.win1 / 10.
        self.assertEqual(win[4,5],0.2)

    def test_radd(self):
        win = 10. + self.win1
        self.assertEqual(win[4,5],12.)

class TestCCD(unittest.TestCase):

    def setUp(self):
        """
        Create two Windows, place into CCDs
        """
        llx1, lly1, xbin, ybin = 12, 22, 2, 3
        dat1  = np.zeros((10,20))
        dat1 += 5
        win1 = ultracam.Window(dat1, llx1, lly1, xbin, ybin) 

        llx2, lly2, xbin, ybin = 120, 220, 2, 3
        dat2  = dat1 + 10.
        win2 = ultracam.Window(dat2, llx2, lly2, xbin, ybin) 

        time = ultracam.Time(56000.2, 5., True, '')

        head = ultracam.Uhead()
        head.add_entry('Object', 'IP Peg', ultracam.ITYPE_STRING,'Name of target')
        self.ccd = ultracam.CCD([win1,win2],time,900,800,True,head)


    def test_index(self):
        self.assertEqual(self.ccd[1][4,5],15.)

    def test_plot(self):
        def ok():
            pg.pgopen('/null')
            self.ccd.plot(5,95)
            pg.pgclos()
            return True
        self.assertTrue(ok())

if __name__ == '__main__':
    suite = unittest.TestLoader().loadTestsFromTestCase(TestWindow)
    unittest.TextTestRunner(verbosity=2).run(suite)

    suite = unittest.TestLoader().loadTestsFromTestCase(TestCCD)
    unittest.TextTestRunner(verbosity=2).run(suite)
