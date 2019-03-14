using WaveletMatrices: WaveletMatrix, getindex, rank, select, freq
using IndexableBitVectors
using Test
using Random

Random.seed!(12345)

function seq_select(a, vec, j)
    if j ≤ 0
        return 0
    end
    for i in 1:length(vec)
        if vec[i] == a
            j -= 1
        end
        if j == 0
            return i
        end
    end
    return 0
end

@testset "WaveletMatrix" begin
    @testset "ordered" begin
        wm = WaveletMatrix([0x00, 0x01, 0x02, 0x03])
        @test length(wm) == 4
        # getindex
        @test wm[1] == 0x00
        @test wm[2] == 0x01
        @test wm[3] == 0x02
        @test wm[4] == 0x03
        @test typeof(wm[4]) == UInt8
        # 0x00
        @test rank(0x00, wm, 0) == 0
        @test rank(0x00, wm, 1) == 1
        @test rank(0x00, wm, 2) == 1
        @test rank(0x00, wm, 3) == 1
        @test rank(0x00, wm, 4) == 1
        @test rank(0x00, wm, 5) == 1
        # 0x01
        @test rank(0x01, wm, 0) == 0
        @test rank(0x01, wm, 1) == 0
        @test rank(0x01, wm, 2) == 1
        @test rank(0x01, wm, 3) == 1
        @test rank(0x01, wm, 4) == 1
        @test rank(0x01, wm, 5) == 1
        # 0x02
        @test rank(0x02, wm, 0) == 0
        @test rank(0x02, wm, 1) == 0
        @test rank(0x02, wm, 2) == 0
        @test rank(0x02, wm, 3) == 1
        @test rank(0x02, wm, 4) == 1
        @test rank(0x02, wm, 5) == 1
        # 0x03
        @test rank(0x03, wm, 0) == 0
        @test rank(0x03, wm, 1) == 0
        @test rank(0x03, wm, 2) == 0
        @test rank(0x03, wm, 3) == 0
        @test rank(0x03, wm, 4) == 1
        @test rank(0x03, wm, 5) == 1

        @test select(0x00, wm, 0) == 0
        @test select(0x01, wm, 0) == 0
        @test select(0x02, wm, 0) == 0
        @test select(0x03, wm, 0) == 0

        @test select(0x00, wm, 1) == 1
        @test select(0x01, wm, 1) == 2
        @test select(0x02, wm, 1) == 3
        @test select(0x03, wm, 1) == 4

        @test select(0x00, wm, 2) == 0
        @test select(0x01, wm, 2) == 0
        @test select(0x02, wm, 2) == 0
        @test select(0x03, wm, 2) == 0

        @test freq(0x00, wm, 1, 4) == 1
        @test freq(0x00, wm, 2, 4) == 0
        @test freq(0x01, wm, 1, 2) == 1
        @test freq(0x01, wm, 1, 1) == 0
        @test freq(0x02, wm, 4, 2) == 0
        @test freq(0x02, wm, 3, 3) == 1
        @test freq(0x03, wm, 4, 1) == 0
        @test freq(0x03, wm, 4, 4) == 1
    end

    @testset "unordered" begin
        wm = WaveletMatrix([0x01, 0x03, 0x02, 0x00])
        @test length(wm) == 4
        # getindex
        @test wm[1] == 0x01
        @test wm[2] == 0x03
        @test wm[3] == 0x02
        @test wm[4] == 0x00
        # 0x00
        @test rank(0x00, wm, 0) == 0
        @test rank(0x00, wm, 1) == 0
        @test rank(0x00, wm, 2) == 0
        @test rank(0x00, wm, 3) == 0
        @test rank(0x00, wm, 4) == 1
        @test rank(0x00, wm, 5) == 1
        # 0x01
        @test rank(0x01, wm, 0) == 0
        @test rank(0x01, wm, 1) == 1
        @test rank(0x01, wm, 2) == 1
        @test rank(0x01, wm, 3) == 1
        @test rank(0x01, wm, 4) == 1
        @test rank(0x01, wm, 5) == 1
        # 0x02
        @test rank(0x02, wm, 0) == 0
        @test rank(0x02, wm, 1) == 0
        @test rank(0x02, wm, 2) == 0
        @test rank(0x02, wm, 3) == 1
        @test rank(0x02, wm, 4) == 1
        @test rank(0x02, wm, 5) == 1
        # 0x03
        @test rank(0x03, wm, 0) == 0
        @test rank(0x03, wm, 1) == 0
        @test rank(0x03, wm, 2) == 1
        @test rank(0x03, wm, 3) == 1
        @test rank(0x03, wm, 4) == 1
        @test rank(0x03, wm, 5) == 1

        @test select(0x00, wm, 1) == 4
        @test select(0x01, wm, 1) == 1
        @test select(0x02, wm, 1) == 3
        @test select(0x03, wm, 1) == 2

        @test select(0x00, wm, 2) == 0
        @test select(0x01, wm, 2) == 0
        @test select(0x02, wm, 2) == 0
        @test select(0x03, wm, 2) == 0
    end

    @testset "homogeneous" begin
        wm = WaveletMatrix([0x01, 0x01, 0x01, 0x01])
        @test length(wm) == 4
        # getindex
        for i in 1:4
            @test wm[i] == 0x01
        end
        @test rank(0x01, wm, 0) == 0
        @test rank(0x01, wm, 1) == 1
        @test rank(0x01, wm, 2) == 2
        @test rank(0x01, wm, 3) == 3
        @test rank(0x01, wm, 4) == 4
        @test rank(0x01, wm, 5) == 4

        @test select(0x01, wm, 1) == 1
        @test select(0x01, wm, 2) == 2
        @test select(0x01, wm, 3) == 3
        @test select(0x01, wm, 4) == 4
        @test select(0x01, wm, 5) == 0

        @test select(0x00, wm, 1) == 0
        @test select(0x02, wm, 1) == 0

        for i in 1:4, j in 1:4
            @test freq(0x00, wm, i, j) == 0
            @test freq(0x01, wm, i, j) == max(0, j - i + 1)
        end
    end

    @testset "Vector{Uint64}" begin
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
        @test length(wm) == 10
        for i in 1:lastindex(x)
            @test wm[i] == x[i]
        end
        @test typeof(wm[1]) == UInt64
        for v in x, i in 1:length(x)
            @test rank(v, wm, i) == count(v′ -> v′ == v, x[1:i])
        end
        for v in x, j in 1:length(x)+1
            @test select(v, wm, j) == seq_select(v, x, j)
        end
    end

    @testset "2-bit encoding" begin
        x = rand(0x00:0x03, 100)
        wm = WaveletMatrix{2}(x)
        @test length(wm) == 100
        for i in 1:lastindex(x)
            @test wm[i] == x[i]
        end
        for a in 0x00:0x03, i in 1:100
            @test rank(a, wm, i) == count(a′ -> a′ == a, x[1:i])
        end
        for a in 0x00:0x03, j in 1:100+1
            @test select(a, wm, j) == seq_select(a, x, j)
        end
    end

    @testset "17-bit encoding" begin
        x = rand(0x00000000:0x00000011, 500)
        wm = WaveletMatrix{17}(x)
        @test length(wm) == 500
        for i in 1:lastindex(x)
            @test wm[i] == x[i]
        end
        for a in 0x00000000:0x00000011, i in 1:500
            @test rank(a, wm, i) == count(a′ -> a′ == a, x[1:i])
        end
        for a in 0x00000000:0x00000011, j in 1:500+1
            @test select(a, wm, j) == seq_select(a, x, j)
        end
    end

    @testset "destructive" begin
        x = rand(0x00:0x03, 1000)
        wm = WaveletMatrix{2}(x, destructive=false)
        @test length(x) == 1000
        wm = WaveletMatrix{2}(x, destructive=true)
        @test length(x) == 0
    end

    @testset "random" begin
        # there is at least one duplicated byte (pigeonhole principle)
        len = typemax(UInt8) + 1
        bytes = rand(UInt8, len)
        wm = WaveletMatrix(copy(bytes))
        for i in 1:len
            @test wm[i] == bytes[i]
        end
        for byte in 0x00:0xff, i in 1:len
            @test rank(byte, wm, i) == count(b -> b == byte, bytes[1:i])
        end
        for byte in 0x00:0xff, j in 0:len+1
            @test select(byte, wm, j) == seq_select(byte, bytes, j)
        end
        for byte in 0x00:0xff
            @test freq(byte, wm, len, 1) == 0
            @test freq(byte, wm, 1, len) == count(b -> b == byte, bytes)
            i, j = rand(1:len), rand(1:len)
            @test freq(byte, wm, i, j) == count(b -> b == byte, bytes[i:j])
        end
    end
end
