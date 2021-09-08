import unittest
from func import *

class TestStringMethods(unittest.TestCase):

  def test_qf(self):
    value = qf(.5, 10, 20, lower=False)
    self.assertAlmostEqual(value, 0.9662639)

  def test_pf(self):
    value = pf(.5, 10, 20, lower=False)
    self.assertAlmostEqual(value, 0.8701604)

if __name__ == '__main__':
  unittest.main()