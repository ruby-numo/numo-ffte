# FFTE interface for Ruby with Numo::NArray

## FFTE
* Package to compute Discrete Fourier Transforms of
  1-, 2- and 3- dimensional sequences of length (2^p)*(3^q)*(5^r),
  developed by Prof. Daisuke Takahashi at University of Tsukuba.
* FFTE site: [FFTE web site](http://www.ffte.jp/)

## Ruby wrapper
* [GitHub site](https://github.com/ruby-numo/ffte)
* [Tentative API Document](http://ruby-numo.github.io/ffte/ref/frames.html)

## Installation
* FORTRAN compiler is required to build FFTE (gfortran, etc).
* Install [Numo::NArray](https://github.com/ruby-numo/narray)
* Install Numo::FFTE:
```shell
$ git clone git://github.com/ruby-numo/ffte.git
$ cd ffte
$ rake build
$ gem install pkg/numo-ffte-*.gem
```

## License
FFTE License:
```
Copyright(C), 2000-2004, 2008-2014, Daisuke Takahashi
You may use, copy, modify this code for any purpose (include
commercial use) and without fee.
You may distribute this ORIGINAL package.
```
