# HCtest -- Higher Criticism Test

Higher Criticism (HC) test statistics of P-values. 

## Example:
```
import numpy as np
from scipy.stats import norm

n = 1000 #number of samples

X = np.random.norm(n)
pvals = norm.sf(X)

hc = HC(pvals)
hv_val, p_th = HC.HCstar(alpha = 0.25)

print("Higher-Criticism statistic = ", hc_val)
```
