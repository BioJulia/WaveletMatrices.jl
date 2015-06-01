# WaveletMatrices

[![Build Status](https://travis-ci.org/bicycle1885/WaveletMatrices.jl.svg?branch=master)](https://travis-ci.org/bicycle1885/WaveletMatrices.jl)

**This package is pre-alpha, unstable, and inefficient.**

An implementation of "The Wavelet Matrix" (Claude and Navarro) <http://www.dcc.uchile.cl/~gnavarro/ps/spire12.4.pdf>.

At the current moment, `rank(a::Uint8, wm::WaveletMatrices, i::Int)` is the only supported operation.
`getindex` and `select` will come soon.


## Requirements

* [IndexedBitVectors.jl](https://github.com/bicycle1885/IndexedBitVectors.jl)


## Demo

This is a full-text search demonstration using the Wavelet matrix and the FM-index algorithm.

Before continuing with this demonstration, you'll need to install the [SuffixArrays.jl](https://github.com/quinnj/SuffixArrays.jl) package to run the Burrows-Wheeler transform (BWT).

Here we'll use the human chromosome 22 (chr22) as a text:

    $ curl -O http://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz
    $ gzip -d chr22.fa.gz

Start up an iteractive session to search DNA fragments; please wait for a few minutes to build an index, then type some DNA queries:

    $ julia fmindex.jl chr22.fa
    INFO: loading a text
    INFO: text length: 51,304,572 characters
    INFO: building a BWT
    INFO: counting characters
    INFO: building a WaveletMatrix
    query> ACGTTG
    'ACGTTG' occurs 1,959 times (9,692 μs).
    query> TATATTATTTAT
    'TATATTATTTAT' occurs 24 times (48 μs).
    query>

As you notice, the session reports the number of occurrences of the pattern and the elapsed time to search.
If you append ';' at the end of a query, the program uses the `Base.search` function to search the query:

    query> ACGTTG;
    'ACGTTG' occurs 1,959 times (130,851 μs).
    query> TATATTATTTAT;
    'TATATTATTTAT' occurs 24 times (60,350 μs).
    query>

Since `Base.search` is a linear search and does not use the index, '..;' search runs far slower:

    query> AAA
    'AAA' occurs 1,016,863 times (16 μs).
    query> AAA;
    'AAA' occurs 1,016,863 times (127,242 μs).
    query>
