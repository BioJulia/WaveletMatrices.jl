module WaveletMatrices

export WaveletMatrix, getindex, rank, freq

import Base: endof, length, sizeof, getindex

using IndexableBitVectors
import IndexableBitVectors: rank

immutable WaveletMatrix{n,B<:AbstractBitVector}
    bits::NTuple{n,B}
    nzeros::NTuple{n,Int}
    len::Int
end

function WaveletMatrix{B<:AbstractBitVector,T<:Unsigned}(::Type{B}, data::Vector{T}, n::Int)
    @assert 1 ≤ n ≤ sizeof(T) * 8
    len = length(data)
    bits = B[]
    nzeros = Int[]
    # TODO: efficient construction
    data = copy(data)
    data′ = Array(T, len)
    for d in 1:n
        # scan d-th bit
        bits′ = B()
        for i in 1:len
            if (data[i] >> (n - d)) & 1 == 1
                # right
                push!(bits′, 1)
            else
                # left
                push!(bits′, 0)
            end
        end
        nzero = rank0(bits′, len)
        l = r = 1
        for i in 1:len
            if bits′[i]
                # right
                data′[nzero+r] = data[i]
                r += 1
            else
                # left
                data′[l] = data[i]
                l += 1
            end
        end
        push!(bits, bits′)
        push!(nzeros, nzero)
        copy!(data, data′)
    end
    return WaveletMatrix{n,B}(tuple(bits...), tuple(nzeros...), len)
end

WaveletMatrix{T<:Unsigned}(src::Vector{T}, n::Int=sizeof(T) * 8) = WaveletMatrix(CompactBitVector, src, n)
WaveletMatrix(str::ByteString) = WaveletMatrix(str.data)

endof(wm::WaveletMatrix)  = wm.len
length(wm::WaveletMatrix) = wm.len

function sizeof{n}(wm::WaveletMatrix{n})
    s = 0
    # bits
    for d in 1:n
        s += sizeof(wm.bits[d])
    end
    # nzeros
    s += n * sizeof(Int)
    # len
    s += sizeof(Int)
    return s
end

function getindex{n}(wm::WaveletMatrix{n}, i::Integer)
    if i < 0 || endof(wm) < i
        throw(BoundsError(i))
    end
    ret = UInt64(0)
    @inbounds for d in 1:n
        bit = wm.bits[d][i]
        if bit
            i = wm.nzeros[d] + rank1(wm.bits[d], i)
        else
            i = rank0(wm.bits[d], i)
        end
        ret |= convert(UInt64, bit) << (n - d)
    end
    return ret
end

function rank{n}(a::Unsigned, wm::WaveletMatrix{n}, i::Int)
    if i < 0 || endof(wm) < i
        throw(BoundsError(i))
    elseif i == 0
        return 0
    end
    sp = 0
    ep = i
    # scan from the most significant bit to the least significant bit
    # e.g.
    #   0b01000101
    #     =======>
    #   d=1 .... 8
    @inbounds for d in 1:n
        bit = (a >> (n - d)) & 1 == 1
        sp = rank(bit, wm.bits[d], sp)
        ep = rank(bit, wm.bits[d], ep)
        if bit
            sp += wm.nzeros[d]
            ep += wm.nzeros[d]
        end
    end
    return ep - sp
end

rank(a::Unsigned, wm::WaveletMatrix, i::Integer) = rank(a, wm, convert(Int, i))

function freq(c::Unsigned, wm::WaveletMatrix, i::Integer, j::Integer)
    return j < i ? 0 : rank(c, wm, j) - rank(c, wm, i - 1)
end

end # module
