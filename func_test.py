import unittest
from func import *

class TestPwrMethods(unittest.TestCase):

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

  def test_pwr_t_test_two_sample(self):
    value = pwr_t_test(100, .1, .05, type='two-sample')
    self.assertAlmostEqual(value, 0.1083718)

  def test_pwr_t_test_one_sample(self):
    value = pwr_t_test(100, .1, .05, type='one-sample')
    self.assertAlmostEqual(value, 0.1676995)

  def test_pwr_t2n_test(self):
    value = pwr_t2n_test(100, 200, .1, .05)
    self.assertAlmostEqual(value, 0.1286476)

  def test_pwr_chisq_test(self):
    value = pwr_chisq_test(.1, 10, 10, .05)
    self.assertAlmostEqual(value, 0.05287121, places=3)

  def test_pwr_f2_test(self):
    value = pwr_f2_test(10, 20, .1, .05)
    self.assertAlmostEqual(value, 0.1257148, places=4)

  def test_pwr_anova_test(self):
    value = pwr_anova_test(2, 3, .1, .05)
    self.assertAlmostEqual(value, 0.05426753)

  def test_pwr_anova_test_find_n(self):
    value = pwr_anova_test(k=2, f=.01, sig_level=.05, power=.5)
    self.assertAlmostEqual(value, 19206.08, places=1)

  def test_pwr_p_test(self):
    value = pwr_p_test(.6, 10, .05)
    self.assertAlmostEqual(value, 0.4751009)

  def test_pwr_2p_test(self):
    value = pwr_2p_test(.6, 10, .05)
    self.assertAlmostEqual(value, 0.2686618)

  def test_pwr_2p2n_test(self):
    value = pwr_2p2n_test(.6, 10, 20, .05)
    self.assertAlmostEqual(value, 0.3408451)

  def test_pwr_norm_test(self):
    value = pwr_norm_test(.6, 10, .05)
    self.assertAlmostEqual(value, 0.4751009)

  def test_pwr_r_test(self):
    value = pwr_r_test(10, .5, .05)
    self.assertAlmostEqual(value, 0.3290749)

  def test_ES_h(self):
    value = ES_h(.1, .2)
    self.assertAlmostEqual(value, -0.2837941)

  def test_ES_w1(self):
    P0 = np.array([.5, .5])
    P1 = np.array([.6, .4])
    value = ES_w1(P0, P1)
    self.assertAlmostEqual(value, 0.2)

  def test_ES_w2(self):
    P = np.array([[.4, .3], [.2, .1]])
    value = ES_w2(P)
    self.assertAlmostEqual(value, 0.08908708)

if __name__ == '__main__':
  unittest.main()