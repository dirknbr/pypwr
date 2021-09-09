import unittest
from func import *

class TestStringMethods(unittest.TestCase):

  def test_qf(self):
    value = qf(.5, 10, 20, lower=False)
    self.assertAlmostEqual(value, 0.9662639)

  def test_pf(self):
    value = pf(.5, 10, 20, lower=False)
    self.assertAlmostEqual(value, 0.8701604)

  def test_qt(self):
    value = qt(.6, 10, lower=False)
    self.assertAlmostEqual(value, -0.2601848)

  def test_pt(self):
    value = pt(.6, 10, lower=False)
    self.assertAlmostEqual(value, 0.2809276)

  def test_qnorm(self):
    value = qnorm(.6, lower=False)
    self.assertAlmostEqual(value, -0.2533471)

  def test_pnorm(self):
    value = pnorm(.6, lower=False)
    self.assertAlmostEqual(value, 0.2742531)

if __name__ == '__main__':
  unittest.main()