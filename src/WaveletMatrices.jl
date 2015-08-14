module WaveletMatrices

export WaveletMatrix, getindex, rank, select, freq

import Base: endof, length, size, sizeof, getindex, select

using IndexableBitVectors
import IndexableBitVectors: rank

include("build.jl")

immutable WaveletMatrix{n,T<:Unsigned,B<:AbstractBitVector} <: AbstractVector{T}
    bits::NTuple{n,B}
    nzeros::NTuple{n,Int}
    sps::Vector{Int}
    function WaveletMatrix(bits)
        @assert 1 ≤ n ≤ sizeof(T) * 8 ≤ 64
        @assert length(bits) == n
        nzeros = Int[]
        for bv in bits
            @assert length(bits[1]) == length(bv)
            push!(nzeros, rank0(bv, length(bv)))
        end
        if n ≤ 16
            # size of lookup table ≤ 512KiB (= sizeof(Int) * 2^16)
            alphabetsize = 2^n
            sps = Vector{Int}(alphabetsize)
            for a in 0:alphabetsize-1
                sps[a+1] = locate_sp(T(a), bits, nzeros)
            end
        else
            sps = Int[]
        end
        return new(tuple(bits...), tuple(nzeros...), sps)
    end
end

const default_bitvector = SucVector

function Base.call{n,T<:Unsigned}(::Type{WaveletMatrix{n}}, data::AbstractVector{T})
    bits = build(default_bitvector, data, n)
    return WaveletMatrix{n,T,default_bitvector}(bits)
end

function Base.call{T<:Unsigned}(::Type{WaveletMatrix}, data::AbstractVector{T}, n::Integer=sizeof(T) * 8)
    bits = build(default_bitvector, data, n)
    return WaveletMatrix{n,T,default_bitvector}(bits)
end

length(wm::WaveletMatrix) = length(wm.bits[1])
size(wm::WaveletMatrix) = (length(wm),)

function sizeof{n}(wm::WaveletMatrix{n})
    s = 0
    for d in 1:n
        s += sizeof(wm.bits[d])
    end
    s += sizeof(wm.nzeros)
    s += sizeof(wm.sps)
    return s
end

@inline function getindex{n,T}(wm::WaveletMatrix{n,T}, i::Int)
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

@inline getindex{n,T}(wm::WaveletMatrix{n,T}, i::Integer) = getindex(wm, convert(Int, i))

function rank{n}(a::Unsigned, wm::WaveletMatrix{n}, i::Int)
    if i < 0 || endof(wm) < i
        throw(BoundsError(i))
    elseif i == 0
        return 0
    end
    if !isempty(wm.sps)
        # use precomputed sp
        sp = wm.sps[a+1]
        ep = i
        @inbounds for d in 1:n
            bit = (a >> (n - d)) & 1 == 1
            ep = rank(bit, wm.bits[d], ep)
            if bit
                ep += wm.nzeros[d]
            end
        end
    else
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
    end
    return ep - sp
end

rank(a::Unsigned, wm::WaveletMatrix, i::Integer) = rank(a, wm, convert(Int, i))

function select(a::Unsigned, wm::WaveletMatrix, j::Int)
    if j ≤ 0
        return 0
    end
    # binary search: j ∈ (rank(l), rank(u)]
    l = 0
    u = length(wm)
    rank_u = rank(a, wm, u)
    while j ≤ rank_u
        m = div(l + u, 2)
        if l == m
            return u
        end
        rank_m = rank(a, wm, m)
        if j ≤ rank_m
            u = m
            rank_u = rank_m
        else
            l = m
        end
    end
    return 0
end

select(a::Unsigned, wm::WaveletMatrix, j::Integer) = select(a, wm, convert(Int, j))

function freq(a::Unsigned, wm::WaveletMatrix, i::Integer, j::Integer)
    return j < i ? 0 : rank(a, wm, j) - rank(a, wm, i - 1)
end

end # module
