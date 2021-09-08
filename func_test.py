import unittest
from func import qf

class TestStringMethods(unittest.TestCase):

  def test_qf(self):
    value = qf(.5, 10, 20, lower=False)
    self.assertAlmostEqual(value, 0.9662639)


if __name__ == '__main__':
  unittest.main()