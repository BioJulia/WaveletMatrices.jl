module WaveletMatrices

export WaveletMatrix, rank

import Base: endof, length, sizeof, getindex

using IndexableBitVectors
import IndexableBitVectors: rank

immutable WaveletMatrix{B<:AbstractBitVector,N}
    bits::NTuple{N,B}
    nzeros::NTuple{N,Int}
    len::Int
end

function WaveletMatrix{B<:AbstractBitVector,T<:Unsigned}(::Type{B}, data::Vector{T}, N::Int)
    @assert 1 ≤ N ≤ sizeof(T) * 8
    len = length(data)
    bits = B[]
    nzeros = Int[]
    # TODO: efficient construction
    data = copy(data)
    data′ = Array(T, len)
    for d in 1:N
        # scan d-th bit
        bits′ = B()
        for i in 1:len
            if (data[i] >> (N - d)) & 1 == 1
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
    return WaveletMatrix{CompactBitVector,N}(tuple(bits...), tuple(nzeros...), len)
end

WaveletMatrix{T}(src::Vector{T}, N::Int=sizeof(T) * 8) = WaveletMatrix(CompactBitVector, src, N)
WaveletMatrix(str::ByteString) = WaveletMatrix(str.data)

endof(wm::WaveletMatrix)  = wm.len
length(wm::WaveletMatrix) = wm.len

function sizeof{B,N}(wm::WaveletMatrix{B,N})
    s = 0
    # bits
    for n in 1:N
        s += sizeof(wm.bits[n])
    end
    # nzeros
    s += N * sizeof(Int)
    # len
    s += sizeof(Int)
    return s
end

function getindex{B,N}(wm::WaveletMatrix{B,N}, i::Integer)
    ret = zero(UInt64)
    for d in 1:N
        bit = wm.bits[d][i]
        if bit
            i = wm.nzeros[d] + rank1(wm.bits[d], i)
        else
            i = rank0(wm.bits[d], i)
        end
        ret |= convert(UInt64, bit) << (N - d)
    end
    return ret
end

function rank{B<:AbstractBitVector,N}(a::Unsigned, wm::WaveletMatrix{B,N}, i::Int)
    @assert 0 ≤ i ≤ endof(wm)
    if i == 0
        return 0
    end
    sp = 0
    ep = i
    # scan from the most significant bit to the least significant bit
    # e.g.
    #   0b01000101
    #     =======>
    #   d=1 .... 8
    @inbounds for d in 1:N
        bit = (a >> (N - d)) & 1 == 1
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
rank(c::Char, wm::WaveletMatrix, i::Integer) = rank(convert(Uint8, c), wm, i)

end # module
