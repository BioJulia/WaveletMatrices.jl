# WaveletMatrices

[![Build Status](https://travis-ci.org/BioJulia/WaveletMatrices.jl.svg?branch=master)](https://travis-ci.org/BioJulia/WaveletMatrices.jl)
[![GitHub latest release](https://img.shields.io/badge/latest-v0.2.0-blue)](https://github.com/BioJulia/WaveletMatrices.jl)
[![GitHub latest release](https://img.shields.io/badge/stable-v0.2.0-green)](https://github.com/BioJulia/WaveletMatrices.jl/releases)
[![Julia versions](https://img.shields.io/badge/julia-1.0+-blue)](https://julialang.org/)

An implementation of "The Wavelet Matrix" (Claude and Navarro) <http://www.dcc.uchile.cl/~gnavarro/ps/spire12.4.pdf>.

The wavelet matrix is a data structure to represent an immutable sequence of
unsigned integers that supports some queries efficiently.

The `WaveletMatrix` type is defined as follows:

```julia
struct WaveletMatrix{w,T<:Unsigned,B<:AbstractBitVector} <: AbstractVector{T}
```

where

* `w`: the number of bits required to encode the unsigned integers (elements)
* `T`: the type of elements
* `B`: the type of bit vectors used to store elements.

To efficiently pack a sequence of unsigned integers, `w` should be as small as possible but enough to encode those integers.
For example, if you want to store integers between 0x00 and 0x03 (i.e. four distinct integers), setting `w = 2 (= ceil(log2(4)))` is the best choice.

The basic operations available on the wavelet matrix are:

* `getindex(wm::WaveletMatrix, i::Integer)`: Return `i`-th element of `wm`.
* `rank(a::Unsigned, wm::WaveletMatrix, i::Integer)`: Count the number of `a`'s occurrences in `wm` between `1:i`.
* `select(a::Unsigned, wm::WaveletMatrix, j::Integer)`: Return the position of the `j`-th occurrence of `a` in `wm`.


## Example

```julia
data = rand(0x00:0x03, 100_000)
wm = WaveletMatrix{2}(data)

wm[129]
rank(0x02, wm, 5555)
partialsort(0x01, wm, 9876)
```
