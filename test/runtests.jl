using WaveletMatrices
using IndexableBitVectors
using Base.Test
using FactCheck

srand(12345)

facts("WaveletMatrix") do
    context("ordered") do
        wm = WaveletMatrix([0x00, 0x01, 0x02, 0x03])
        @fact length(wm) --> 4
        # getindex
        @fact wm[1] --> 0x00
        @fact wm[2] --> 0x01
        @fact wm[3] --> 0x02
        @fact wm[4] --> 0x03
        # 0x00
        @fact rank(0x00, wm, 1) --> 1
        @fact rank(0x00, wm, 2) --> 1
        @fact rank(0x00, wm, 3) --> 1
        @fact rank(0x00, wm, 4) --> 1
        # 0x01
        @fact rank(0x01, wm, 1) --> 0
        @fact rank(0x01, wm, 2) --> 1
        @fact rank(0x01, wm, 3) --> 1
        @fact rank(0x01, wm, 4) --> 1
        # 0x02
        @fact rank(0x02, wm, 1) --> 0
        @fact rank(0x02, wm, 2) --> 0
        @fact rank(0x02, wm, 3) --> 1
        @fact rank(0x02, wm, 4) --> 1
        # 0x03
        @fact rank(0x03, wm, 1) --> 0
        @fact rank(0x03, wm, 2) --> 0
        @fact rank(0x03, wm, 3) --> 0
        @fact rank(0x03, wm, 4) --> 1

        @fact freq(0x00, wm, 1, 4) --> 1
        @fact freq(0x00, wm, 2, 4) --> 0
        @fact freq(0x01, wm, 1, 2) --> 1
        @fact freq(0x01, wm, 1, 1) --> 0
        @fact freq(0x02, wm, 4, 2) --> 0
        @fact freq(0x02, wm, 3, 3) --> 1
        @fact freq(0x03, wm, 4, 1) --> 0
        @fact freq(0x03, wm, 4, 4) --> 1
    end

    context("unordered") do
        wm = WaveletMatrix([0x01, 0x03, 0x02, 0x00])
        @fact length(wm) --> 4
        # getindex
        @fact wm[1] --> 0x01
        @fact wm[2] --> 0x03
        @fact wm[3] --> 0x02
        @fact wm[4] --> 0x00
        # 0x00
        @fact rank(0x00, wm, 1) --> 0
        @fact rank(0x00, wm, 2) --> 0
        @fact rank(0x00, wm, 3) --> 0
        @fact rank(0x00, wm, 4) --> 1
        # 0x01
        @fact rank(0x01, wm, 1) --> 1
        @fact rank(0x01, wm, 2) --> 1
        @fact rank(0x01, wm, 3) --> 1
        @fact rank(0x01, wm, 4) --> 1
        # 0x02
        @fact rank(0x02, wm, 1) --> 0
        @fact rank(0x02, wm, 2) --> 0
        @fact rank(0x02, wm, 3) --> 1
        @fact rank(0x02, wm, 4) --> 1
        # 0x03
        @fact rank(0x03, wm, 1) --> 0
        @fact rank(0x03, wm, 2) --> 1
        @fact rank(0x03, wm, 3) --> 1
        @fact rank(0x03, wm, 4) --> 1
    end

    context("homogeneous") do
        wm = WaveletMatrix([0x01, 0x01, 0x01, 0x01])
        @fact length(wm) --> 4
        # getindex
        for i in 1:4
            @fact wm[i] --> 0x01
        end
        @fact rank(0x01, wm, 1) --> 1
        @fact rank(0x01, wm, 2) --> 2
        @fact rank(0x01, wm, 3) --> 3
        @fact rank(0x01, wm, 4) --> 4

        for i in 1:4, j in 1:4
            @fact freq(0x00, wm, i, j) --> 0
            @fact freq(0x01, wm, i, j) --> max(0, j - i + 1)
        end
    end

    context("Vector{Uint64}") do
        x = [
            0x0000000000000001,
            0xffffffffffffffff,
            0x0000000000000000,
            0x0000000000000001,
            0xffffffffffffffff,
            0xffffff171f410d54,
            0x0000000000000009,
            0x0000000000000000,
            0x1fffff171f410d54,
            0x0000000000000008,
        ]
        wm = WaveletMatrix(x)
        @fact length(wm) --> 10
        for i in 1:endof(x)
            @fact wm[i] --> x[i]
        end
        for v in x, i in 1:length(x)
            @fact rank(v, wm, i) --> count(v′ -> v′ == v, x[1:i])
        end
    end

    context("2-bit encoding") do
        x = rand(0x00:0x03, 100)
        wm = WaveletMatrix{2}(x)
        @fact length(wm) --> 100
        for i in 1:endof(x)
            @fact wm[i] --> x[i]
        end
        for a in 0x00:0x03, i in 1:100
            @fact rank(a, wm, i) --> count(a′ -> a′ == a, x[1:i])
        end
    end

    context("17-bit encoding") do
        x = rand(0x00000000:0x00000011, 500)
        wm = WaveletMatrix{17}(x)
        @fact length(wm) --> 500
        for i in 1:endof(x)
            @fact wm[i] --> x[i]
        end
        for a in 0x00000000:0x00000011, i in 1:500
            @fact rank(a, wm, i) --> count(a′ -> a′ == a, x[1:i])
        end
    end

    context("random") do
        # there is at least one duplicated byte (pigeonhole principle)
        len = typemax(Uint8) + 1
        bytes = rand(Uint8, len)
        wm = WaveletMatrix(copy(bytes))
        for i in 1:len
            @fact wm[i] --> bytes[i]
        end
        for byte in 0x00:0xff, i in 1:len
            @fact rank(byte, wm, i) --> count(b -> b == byte, bytes[1:i])
        end
        for byte in 0x00:0xff
            @fact freq(byte, wm, len, 1) --> 0
            @fact freq(byte, wm, 1, len) --> count(b -> b == byte, bytes)
            i, j = rand(1:len), rand(1:len)
            @fact freq(byte, wm, i, j) --> count(b -> b == byte, bytes[i:j])
        end
    end
end
