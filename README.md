# FFTE Ruby wrapper with Numo::NArray

## FFTE
* is a package to compute Discrete Fourier Transforms of
  1-, 2- and 3- dimensional sequences of length (2^p)*(3^q)*(5^r),
  developed by Prof. Daisuke Takahashi at University of Tsukuba.
* Visit [FFTE web site](http://www.ffte.jp/)

## Ruby wrapper
* [GitHub site](https://github.com/masa16/numo-ffte)
* [Tentative API Document](http://masa16.github.io/ffte/ref/frames.html)

## Installation
* FORTRAN compiler is required to build FFTE (gfortran, etc).
* Install [Numo::NArray](https://github.com/masa16/numo-narray)
* Install Numo::FFTE:
```shell
$ git clone git://github.com/masa16/numo-ffte.git
$ cd numo-ffte
$ rake build
$ gem install pkg/numo-ffte-*.gem
```
