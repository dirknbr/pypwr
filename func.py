
# we keep the general synatx of the functions from R
# dots become hyphens or underscores

import numpy as np
import scipy.stats as stats
from scipy.optimize import root
from statsmodels.stats.power import tt_ind_solve_power, tt_solve_power, FTestPower, GofChisquarePower

def pwr_t_test(n=None, d=None, sig_level=.05, power=None,
  type='two-sample', alternative='two-sided'):
  # greater = larger
  # smaller = less
  if type == 'two-sample':
    return tt_ind_solve_power(effect_size=d, nobs1=n, alpha=sig_level, power=power,
                              alternative=alternative)
  elif type == 'one-sample':
    # https://www.statsmodels.org/stable/generated/statsmodels.stats.power.tt_solve_power.html#statsmodels.stats.power.tt_solve_power
    return tt_solve_power(effect_size=d, nobs=n, alpha=sig_level, power=power, 
                          alternative=alternative)

def pwr_chisq_test(w=None, N=None, df=None, sig_level=.05, power=None):
  # https://www.statsmodels.org/stable/generated/statsmodels.stats.power.GofChisquarePower.solve_power.html#statsmodels.stats.power.GofChisquarePower.solve_power
  test = GofChisquarePower()
  return test.solve_power(effect_size=w, nobs=N, alpha=sig_level,
                          power=power, n_bins=df)

def qf(p, df1, df2, ncp=None, lower=True):
  if not lower:
    p = 1 - p
  if ncp:
    # https://docs.scipy.org/doc/scipy/reference/reference/generated/scipy.stats.ncf.html#scipy.stats.ncf
    return stats.ncf.ppf(p, df1, df2, ncp)
  else:
    return stats.f.ppf(p, df1, df2)

def pf(q, df1, df2, ncp=None, lower=True):
  if ncp:
    value = stats.ncf.cdf(q, df1, df2, ncp)
  else:
    value = stats.f.cdf(q, df1, df2)
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
  _lambda = f2 * (u + v + 1)
  # TODO: fix, returns wrong power
  return qf(pf(sig_level, u, v, lower=False), u, v, _lambda, lower=False)

def pwr_anova_test(k=None, n=None, f=None, sig_level=.05, power=None):
  # TODO: also allow solve for k, n, f, sig_level
  _lambda = k * n * f ** 2
  return pf(qf(sig_level, k - 1, (n - 1) * k, lower=False), k - 1, (n - 1) * k, _lambda, lower=False)

def pwr_p_test(h=None, n=None, sig_level=.05, power=None, alternative='two-sided'):
  # TODO: also allow solve for h, n and sig_level
  if alternative == 'two-sided':
    value = (pnorm(qnorm(sig_level / 2, lower=False) - h * np.sqrt(n), lower=False) +
      pnorm(qnorm(sig_level / 2, lower=True) - h * np.sqrt(n), lower=True))
  elif alternative == 'less':
    value = pnorm(qnorm(sig_level, lower=True) - h * np.sqrt(n), lower=True)
  elif alternative == 'greater':
    value = pnorm(qnorm(sig_level, lower=False) - h * np.sqrt(n), lower=False)    
  return value
  
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
print('qf', qf(.5, 10, 20, lower=False))
print('pf', pf(.5, 10, 20, lower=False))
print('qnorm', qnorm(.6, lower=False))
print('pnorm', pnorm(.6, lower=False))

print('t-test 2 sample', pwr_t_test(100, .1, .05, type='two-sample'))
print('t-test 1 sample', pwr_t_test(100, .1, .05, type='one-sample'))
print('chisq', pwr_chisq_test(.1, 10, 10, .05))
print('f2', pwr_f2_test(10, 20, .1, .05))
print('anova', pwr_anova_test(2, 3, .1, .05))
print('p test', pwr_p_test(.6, 10, .05))

print('ES h', ES_h(.1, .2))
P0 = np.array([.5, .5])
P1 = np.array([.6, .4])
print('ES h1', ES_w1(P0, P1))
P = np.array([[.4, .3], [.2, .1]])
print('ES w2', ES_w2(P))
