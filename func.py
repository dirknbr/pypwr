
# we keep the general synatx of the functions from R
# dots become hyphens or underscores

import numpy as np
import scipy.stats as stats
from scipy.optimize import root, brentq
from statsmodels.stats.power import FTestPower, GofChisquarePower

NRANGE = (2, 1e+09)
SIGRANGE = (1e-10, .9999)

def pwr_t_test(n=None, d=None, sig_level=.05, power=None,
               type='two-sample', alternative='two-sided'):
  tside = 2 if alternative == 'two-sided' else 1
  tsample = 2 if type == 'two-sample' else 1
  ncp = np.sqrt(n / tsample) * d
  nu = (n - 1) * tsample
  if alternative == 'two-sided':
    def _power(n, d, sig_level):
      qu = qt(sig_level / tside, nu, lower=False)
      return pt(qu, nu, ncp=ncp, lower=False) + pt(-qu, nu, ncp=ncp, lower=True)
    drange = (1e-07, 10)
  elif alternative == 'less':
    def _power(n, d, sig_level):
      return pt(qt(sig_level / tside, nu, lower=True), nu, ncp=ncp, lower=True)
    drange = (-10, 5)
  elif alternative == 'greater':
    def _power(n, d, sig_level):
      return pt(qt(sig_level / tside, nu, lower=False), nu, ncp=ncp, lower=False)
    drange = (-5, 10)
  if power is None:
    return _power(n, d, sig_level)
  elif n is None:
    def __power(n):
      return _power(n, d, sig_level) - power
    return brentq(__power, 2, 1e+09)
  elif d is None:
    def __power(d):
      return _power(n, d, sig_level) - power
    return brentq(__power, drange[0], drange[1])
  elif sig_level is None:
    def __power(sig_level):
      return _power(n, d, sig_level) - power
    return brentq(__power, SIGRANGE[0], SIGRANGE[1])

def pwr_t2n_test(n1=None, n2=None, d=None, sig_level=.05, power=None, 
                 alternative='two-sided'):
  tsample = 2
  tside = 2 if alternative == 'two-sided' else 1
  ncp = d * (1 / np.sqrt(1 / n1 + 1 / n2))
  nu = n1 + n2 - 2
  if alternative == 'two-sided':
    def _power(n1, n2, d, sig_level):
      qu = qt(sig_level / tside, nu, lower=False)
      return pt(qu, nu, ncp, lower=False) + pt(-qu, nu, ncp, lower=True)
    drange = (1e-07, 10)
  elif alternative == 'less':
    def _power(n1, n2, d, sig_level):
      return pt(qt(sig_level / tside, nu, lower=True), nu, ncp, lower=True)
    drange = (-10, 5)
  elif alternative == 'greater':
    def _power(n1, n2, d, sig_level):
      return pt(qt(sig_level / tside, nu, lower=False), nu, ncp, lower=False)
    drange = (-5, 10)
  if power is None:
    return _power(n1, n2, d, sig_level)
  elif n1 is None:
    def __power(n1):
      return _power(n1, n2, d, sig_level) - power
    return brentq(__power, NRANGE[0], NRANGE[1])
  elif n2 is None:
    def __power(n2):
      return _power(n1, n2, d, sig_level) - power
    return brentq(__power, NRANGE[0], NRANGE[1])
  elif d is None:
    def __power(d):
      return _power(n1, n2, d, sig_level) - power
    return brentq(__power, drange[0], drange[1])
  elif sig_level is None:
    def __power(sig_level):
      return _power(n1, n2, d, sig_level)
    return brentq(__power, SIGRANGE[0], SIGRANGE[1])

def pwr_chisq_test(w=None, N=None, df=None, sig_level=.05, power=None):
  # https://www.statsmodels.org/stable/generated/statsmodels.stats.power.GofChisquarePower.solve_power.html#statsmodels.stats.power.GofChisquarePower.solve_power
  test = GofChisquarePower()
  return test.solve_power(effect_size=w, nobs=N, alpha=sig_level,
                          power=power, n_bins=df)

def qf(p, df1, df2, ncp=None, lower=True):
  if not lower:
    p = 1 - p
  if ncp is None: ncp = 0
  return stats.ncf.ppf(p, df1, df2, ncp)

def pf(q, df1, df2, ncp=None, lower=True):
  if ncp is None: ncp = 0
  value = stats.ncf.cdf(q, df1, df2, ncp)
  if not lower:
    return 1 - value
  else:
    return value

def qt(p, df, ncp=None, lower=True):
  if not lower:
    p = 1 - p
  if ncp is None: ncp = 0
  return stats.nct.ppf(p, df, ncp)

def pt(q, df, ncp=None, lower=True):
  if ncp is None: ncp = 0
  value = stats.nct.cdf(q, df, ncp)
  if not lower:
    return 1 - value
  else:
    return value

def qnorm(p, lower=True):
  if not lower:
    p = 1 - p
  return stats.norm.ppf(p)

def pnorm(q, lower=True):
  if not lower:
    return 1 - stats.norm.cdf(q)
  else:
    return stats.norm.cdf(q)
  
def pwr_f2_test(u=None, v=None, f2=None, sig_level=.05, power=None):
  def _power(u, v, f2, sig_level):
    _lambda = f2 * (u + v + 1)
    return pf(qf(sig_level, u, v, lower=False), u, v, _lambda, lower=False)
  if power is None:
    return _power(u, v, f2, sig_level)
  elif u is None:
    def __power(u):
      return _power(u, v, f2, sig_level) - power
    return brentq(__power, 1, 100)
  elif v is None:
    def __power(v):
      return _power(u, v, f2, sig_level) - power
    return brentq(__power, 1, 1e+09)
  elif f2 is None:
    def __power(f2):
      return _power(u, v, f2, sig_level) - power
    return brentq(__power, 1e-07, 1e+07)
  elif sig_level is None:
    def __power(sig_level):
      return _power(u, v, f2, sig_level) - power
    return brentq(__power, SIGRANGE[0], SIGRANGE[1])

def pwr_anova_test(k=None, n=None, f=None, sig_level=.05, power=None):
  def _power(k, n, f, sig_level):
    _lambda = k * n * f ** 2
    return pf(qf(sig_level, k - 1, (n - 1) * k, lower=False), k - 1, (n - 1) * k, _lambda, lower=False)
  if power is None:
    return _power(k, n, f, sig_level)
  elif k is None:
    def __power(k):
      return _power(k, n, f, sig_level) - power
    return brentq(__power, 2, 100)
  elif n is None:
    def __power(n):
      return _power(k, n, f, sig_level) - power
    return brentq(__power, NRANGE[0], NRANGE[1])
  elif f is None:
    def __power(f):
      return _power(k, n, f, sig_level) - power
    return brentq(__power, 1e-07, 1e+07)
  elif sig_level is None:
    def __power(sig_level):
      return _power(k, n, f, sig_level) - power
    return brentq(__power, SIGRANGE[0], SIGRANGE[1])

def pwr_p_test(h=None, n=None, sig_level=.05, power=None, alternative='two-sided'):
  if alternative == 'two-sided':
    def _power(h, n, sig_level):
      return (pnorm(qnorm(sig_level / 2, lower=False) - h * np.sqrt(n), lower=False) +
              pnorm(qnorm(sig_level / 2, lower=True) - h * np.sqrt(n), lower=True))
    hrange = (1e-10, 10)
  elif alternative == 'less':
    def _power(h, n, sig_level):
      return pnorm(qnorm(sig_level, lower=True) - h * np.sqrt(n), lower=True)
    hrange = (-10, 5)
  elif alternative == 'greater':
    def _power(h, n, sig_level):
      return pnorm(qnorm(sig_level, lower=False) - h * np.sqrt(n), lower=False)    
    hrange = (-5, 10)
  if power is None:
    return _power(h, n, sig_level)
  elif h is None:
    def __power(h):
      return _power(h, n, sig_level) - power
    return brentq(__power, hrange[0], hrange[1])
  elif n is None:
    def __power(n):
      return _power(h, n, sig_level) - power
    return brentq(__power, NRANGE[0], NRANGE[1])
  elif sig_level is None:
    def __power(sig_level):
      return _power(h, n, sig_level) - power
    return brentq(__power, SIGRANGE[0], SIGRANGE[1])


def pwr_2p_test(h=None, n=None, sig_level=.05, power=None, alternative='two-sided'):
  return pwr_p_test(h=h, n=n / 2, sig_level=sig_level, power=power, alternative=alternative)

def pwr_2p2n_test(h=None, n1=None, n2=None, sig_level=.05, power=None, alternative='two-sided'):
  n = 2 * (n1 * n2) / (n1 + n2)
  return pwr_2p_test(h=h, n=n, sig_level=sig_level, power=power, alternative=alternative)

def pwr_norm_test(d=None, n=None, sig_level=.05, power=None, alternative='two-sided'):
  return pwr_p_test(d, n, sig_level, power, alternative)

def pwr_r_test(n=None, r=None, sig_level=.05, power=None, alternative='two-sided'):
  if alternative == 'two-sided':
    def _power(n, r, sig_level):
      ttt = qt(sig_level / 2, df=n - 2, lower=False)
      rc = np.sqrt(ttt ** 2 / (ttt ** 2 + n - 2))
      zr = np.arctanh(r) + r / (2 * (n - 1))
      zrc = np.arctanh(rc)
      return pnorm((zr - zrc) * np.sqrt(n - 3)) + pnorm((-zr - zrc) * np.sqrt(n - 3))
  elif alternative == 'less':
    def _power(n, r, sig_level):
      r = -r
      ttt = qt(sig_level, df=n - 2, lower=False)
      rc = np.sqrt(ttt ** 2 / (ttt ** 2 + n - 2))
      zr = np.arctanh(r) + r / (2 * (n - 1))
      zrc = np.arctanh(rc)
      return pnorm((zr - zrc) * np.sqrt(n - 3))
  elif alternative == 'greater':
    def _power(n, r, sig_level):
      ttt = qt(sig_level, df=n - 2, lower=False)
      rc = np.sqrt(ttt ** 2 / (ttt ** 2 + n - 2))
      zr = np.arctanh(r) + r / (2 * (n - 1))
      zrc = np.arctanh(rc)
      return pnorm((zr - zrc) * np.sqrt(n - 3))
  if power is None:
    return _power(n, r, sig_level)
  elif n is None:
    def __power(n):
      return _power(n, r, sig_level) - power
    return brentq(__power, 4, 1e+09)
  elif r is None:
    def __power(r):
      return _power(n, r, sig_level) - power
    if alternative == 'two-sided':
      return brentq(__power, SIGRANGE[0], SIGRANGE[1])
    else:
      return brentq(__power, -.9999, .9999)
  elif sig_level is None:
    def __power(sig_level):
      return _power(n, r, sig_level) - power
    return brentq(__power, SIGRANGE[0], SIGRANGE[1])


def ES_h(p1, p2):
  return 2 * np.arcsin(np.sqrt(p1)) - 2 * np.arcsin(np.sqrt(p2))

def ES_w1(P0, P1):
  assert len(P0) == len(P1)
  return np.sqrt(sum((P1 - P0) ** 2 / P0))

def ES_w2(P):
  assert len(P.shape) == 2 and np.isclose(P.sum(), 1)
  pi = P.sum(axis=1)
  pj = P.sum(axis=0)
  P0 = pi.reshape((len(pi), 1)).dot(pj.reshape((1, len(pj))))
  return np.sqrt(np.sum((P - P0) ** 2 / P0))

# test it
# print('qf', qf(.5, 10, 20, lower=False))
# print('pf', pf(.5, 10, 20, lower=False))
# print('qt', qt(.6, 10, lower=False))
# print('pt', pt(.6, 10, lower=False))
# print('qnorm', qnorm(.6, lower=False))
# print('pnorm', pnorm(.6, lower=False))

# print('t-test 2 sample', pwr_t_test(100, .1, .05, type='two-sample'))
# print('t-test 1 sample', pwr_t_test(100, .1, .05, type='one-sample'))
# print('t-test 2n', pwr_t2n_test(100, 200, .1, .05))
# print('chisq', pwr_chisq_test(.1, 10, 10, .05))
# print('f2', pwr_f2_test(10, 20, .1, .05))

# print('anova', pwr_anova_test(2, 3, .1, .05))
# print('anova n', pwr_anova_test(k=2, n=None, f=.01, sig_level=.05, power=.5))

# print('p test', pwr_p_test(.6, 10, .05))
# print('2p test', pwr_2p_test(.6, 10, .05))
# print('2p2n test', pwr_2p2n_test(.6, 10, 20, .05))
# print('norm test', pwr_norm_test(.6, 10, .05))
# print('r test', pwr_r_test(10, .5, .05))

# print('ES h', ES_h(.1, .2))
# P0 = np.array([.5, .5])
# P1 = np.array([.6, .4])
# print('ES h1', ES_w1(P0, P1))
# P = np.array([[.4, .3], [.2, .1]])
# print('ES w2', ES_w2(P))
