# TracyWidom

[![Build Status](https://travis-ci.org/quit()/TracyWidom.jl.svg?branch=master)](https://travis-ci.org/quit()/TracyWidom.jl)

[![Coverage Status](https://coveralls.io/repos/quit()/TracyWidom.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/quit()/TracyWidom.jl?branch=master)

[![codecov.io](http://codecov.io/github/quit()/TracyWidom.jl/coverage.svg?branch=master)](http://codecov.io/github/quit()/TracyWidom.jl?branch=master)

This is a Julia implementation of Bornemann's method of evaluating the Fredholm determinant associated with the Tracy-Widom distribution.

Provides two functions `F1` and `F2` giving the Tracy-Widom CDFs for GOE and GUE respectively. The keyword argument `num_points` specifies the number of points in the Gauss-Legendre quadrature, with a default of 25

# References
- Folkmar Bornemann
    "On the numerical evaluation of Fredholm determinants",
    *Mathematics of Computation*,
    Volume 79, Number 270, April 2010, Pages 871–915
  [[pdf]](https://www.ams.org/journals/mcom/2010-79-270/S0025-5718-09-02280-7/S0025-5718-09-02280-7.pdf)
  [[doi]](https://doi.org/10.1090/S0025-5718-09-02280-7 )
