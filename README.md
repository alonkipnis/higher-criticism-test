# HCtest -- Higher Criticism Test

Higher Criticism (HC) test statistics for testing the global significance of many independent hypotheses. As an input, 
the test receives a list of P-values and returns the HC test statistics. See (Donoho & Jin 2004)

## Example:
```
from scipy.stats import norm

n = 1000 #number of samples

X = norm.rvs(size=n)
pvals = norm.sf(X)

hc = HC(pvals)
hc_val, p_th = hc.HCstar(gamma = 0.25)

print("Higher-Criticism statistic = ", hc_val)
```
