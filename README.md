# TracyWidom

[![Build Status](https://travis-ci.org/AustenLamacraft/TracyWidom.jl.svg?branch=master)](https://travis-ci.org/AustenLamacraft/TracyWidom.jl)

[![Coverage Status](https://coveralls.io/repos/AustenLamacraft/TracyWidom.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/AustenLamacraft/TracyWidom.jl?branch=master)

[![codecov.io](http://codecov.io/github/AustenLamacraft/TracyWidom.jl/coverage.svg?branch=master)](http://codecov.io/github/AustenLamacraft/TracyWidom.jl?branch=master)

This is a Julia implementation of Bornemann's method of evaluating the Fredholm determinant associated with the Tracy-Widom distribution.

Provides two functions `F1` and `F2` giving the Tracy-Widom CDFs for GOE and GUE respectively. The keyword argument `num_points` specifies the number of points in the Gauss-Legendre quadrature, with a default of 25.

## Sample usage

```
using TracyWidom
F1(0)
0.8319080662029523
F2(0)
0.9693728283552623
F1([0.1,0.2])
2-element Array{Float64,1}:
 0.849358
 0.865444
```

## References
- Folkmar Bornemann
    "On the numerical evaluation of Fredholm determinants",
    *Mathematics of Computation*,
    Volume 79, Number 270, April 2010, Pages 871–915
  [[pdf]](https://www.ams.org/journals/mcom/2010-79-270/S0025-5718-09-02280-7/S0025-5718-09-02280-7.pdf)
  [[doi]](https://doi.org/10.1090/S0025-5718-09-02280-7 )
