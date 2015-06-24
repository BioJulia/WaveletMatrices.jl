using WaveletMatrices
using IndexableBitVectors
using Base.Test
using FactCheck

srand(12345)

facts("WaveletMatrix") do
    context("ordered") do
        wm = WaveletMatrix([0x00, 0x01, 0x02, 0x03])
        # 0x00
        @fact rank(0x00, wm, 1) => 1
        @fact rank(0x00, wm, 2) => 1
        @fact rank(0x00, wm, 3) => 1
        @fact rank(0x00, wm, 4) => 1
        # 0x01
        @fact rank(0x01, wm, 1) => 0
        @fact rank(0x01, wm, 2) => 1
        @fact rank(0x01, wm, 3) => 1
        @fact rank(0x01, wm, 4) => 1
        # 0x02
        @fact rank(0x02, wm, 1) => 0
        @fact rank(0x02, wm, 2) => 0
        @fact rank(0x02, wm, 3) => 1
        @fact rank(0x02, wm, 4) => 1
        # 0x03
        @fact rank(0x03, wm, 1) => 0
        @fact rank(0x03, wm, 2) => 0
        @fact rank(0x03, wm, 3) => 0
        @fact rank(0x03, wm, 4) => 1
    end

    context("unordered") do
        wm = WaveletMatrix([0x01, 0x03, 0x02, 0x00])
        # 0x00
        @fact rank(0x00, wm, 1) => 0
        @fact rank(0x00, wm, 2) => 0
        @fact rank(0x00, wm, 3) => 0
        @fact rank(0x00, wm, 4) => 1
        # 0x01
        @fact rank(0x01, wm, 1) => 1
        @fact rank(0x01, wm, 2) => 1
        @fact rank(0x01, wm, 3) => 1
        @fact rank(0x01, wm, 4) => 1
        # 0x02
        @fact rank(0x02, wm, 1) => 0
        @fact rank(0x02, wm, 2) => 0
        @fact rank(0x02, wm, 3) => 1
        @fact rank(0x02, wm, 4) => 1
        # 0x03
        @fact rank(0x03, wm, 1) => 0
        @fact rank(0x03, wm, 2) => 1
        @fact rank(0x03, wm, 3) => 1
        @fact rank(0x03, wm, 4) => 1
    end

    context("homogeneous") do
        wm = WaveletMatrix([0x01, 0x01, 0x01, 0x01])
        @fact rank(0x01, wm, 1) => 1
        @fact rank(0x01, wm, 2) => 2
        @fact rank(0x01, wm, 3) => 3
        @fact rank(0x01, wm, 4) => 4
    end

    context("ASCIIString") do
        s = "abracadabra"
        wm = WaveletMatrix(s)
        for char in 'a':'z', i in 1:length(s)
            @fact rank(char, wm, i) => count(c -> c == char, s[1:i])
        end
    end

    context("Vector{Uint64}") do
        x = convert(
            Vector{Uint64},
            [1, -1000203023020, 53929202802, 9, 0, typemin(Int64), 2, 2, typemax(Int64), -5, 9]
        )
        wm = WaveletMatrix(x)
        for v in x, i in 1:length(x)
            @fact rank(v, wm, i) => count(v′ -> v′ == v, x[1:i])
        end
    end

    context("2-bit encoding") do
        x = rand(0x00:0x03, 100)
        wm = WaveletMatrix(x, 2)
        for a in 0x00:0x03, i in 1:100
            @fact rank(a, wm, i) => count(a′ -> a′ == a, x[1:i])
        end
    end

    context("17-bit encoding") do
        x = rand(0x00000000:0x00000011, 500)
        wm = WaveletMatrix(x, 17)
        for a in 0x00000000:0x00000011, i in 1:500
            @fact rank(a, wm, i) => count(a′ -> a′ == a, x[1:i])
        end
    end

    context("random") do
        # there is at least one duplicated byte (pigeonhole principle)
        len = typemax(Uint8) + 1
        bytes = rand(Uint8, len)
        wm = WaveletMatrix(copy(bytes))
        for byte in 0x00:0xff, i in 1:len
            @fact rank(byte, wm, i) => count(b -> b == byte, bytes[1:i])
        end
    end
end
