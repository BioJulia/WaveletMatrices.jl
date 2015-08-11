module WaveletMatrices

export WaveletMatrix, getindex, rank, freq

import Base: endof, length, size, sizeof, getindex

using IndexableBitVectors
import IndexableBitVectors: rank

include("build.jl")

immutable WaveletMatrix{n,T<:Unsigned,B<:AbstractBitVector} <: AbstractVector{T}
    bits::NTuple{n,B}
    nzeros::NTuple{n,Int}
    function WaveletMatrix(data::Vector{T})
        @assert 1 ≤ n ≤ sizeof(T) * 8
        bits, nzeros = build(B, data, n)
        new(bits, nzeros)
    end
end

function Base.call{n,T<:Unsigned}(::Type{WaveletMatrix{n}}, data::Vector{T})
    return WaveletMatrix{n,T,CompactBitVector}(data)
end

function Base.call{T<:Unsigned}(::Type{WaveletMatrix}, data::Vector{T}, n::Integer=sizeof(T) * 8)
    return WaveletMatrix{n}(data)
end

length(wm::WaveletMatrix) = length(wm.bits[1])
size(wm::WaveletMatrix) = (length(wm),)

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

function getindex{n,T}(wm::WaveletMatrix{n,T}, i::Integer)
    if i < 0 || endof(wm) < i
        throw(BoundsError(i))
    end
    ret = T(0)
    @inbounds for d in 1:n
        bits = wm.bits[d]
        bit = bits[i]
        ret = ret << 1 | bit
        if d == n
            return ret
        end
        if bit
            i = wm.nzeros[d] + rank1(bits, i)
        else
            i = rank0(bits, i)
        end
    end
end

function rank{n}(a::Unsigned, wm::WaveletMatrix{n}, i::Int)
    if i < 0 || endof(wm) < i
        throw(BoundsError(i))
    elseif i == 0
        return 0
    end
    sp, ep = 0, i
    # scan from the most significant bit to the least significant bit
    @inbounds for d in 1:n
        bit = (a >> (n - d)) & 1 == 1
        sp = rank(bit, wm.bits[d], sp)
        ep = rank(bit, wm.bits[d], ep)
        if bit
            nz = wm.nzeros[d]
            sp += nz
            ep += nz
        end
    end
    return ep - sp
end

rank(a::Unsigned, wm::WaveletMatrix, i::Integer) = rank(a, wm, convert(Int, i))

function freq(c::Unsigned, wm::WaveletMatrix, i::Integer, j::Integer)
    return j < i ? 0 : rank(c, wm, j) - rank(c, wm, i - 1)
end

end # module
