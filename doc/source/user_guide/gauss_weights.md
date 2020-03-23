gauss_weights.py
================

 - Computes the Gaussian weights as a function of degree  
 - A normalized version of [Christopher Jekeli's Gaussian averaging function](http://www.geology.osu.edu/~jekeli.1/OSUReports/reports/report_327.pdf)  

#### Calling Sequence
```
from gravity_toolkit.gauss_weights import gauss_weights
wl = 2.0*np.pi*gauss_weights(hw,LMAX)
```

#### Inputs
 1. `hw`: Gaussian smoothing radius in km  
 2. `LMAX`: Upper bound of Spherical Harmonic Degrees  

#### Outputs
 - `wl`: Gaussian weights for each degree `l`
