tsamplitude.py
==============

- Calculate the amplitude and phase of a harmonic function from calculated sine and cosine of a series of measurements

#### Calling Sequence
```python
from gravity_toolkit.tsamplitude import tsamplitude
ampl,ph = tsamplitude(bsin,bcos)
```
[Source code](https://github.com/tsutterley/read-GRACE-harmonics/blob/main/gravity_toolkit/tsamplitude.py)

#### Inputs
- `bsin`: amplitude of the calculated sine values
- `bcos`: amplitude of the calculated cosine values

#### Outputs
- `ampl`: amplitude from the harmonic functionss
- `ph`: phase from the harmonic functions in degrees
