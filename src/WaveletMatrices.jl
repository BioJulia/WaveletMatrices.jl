module WaveletMatrices

export WaveletMatrix, rank

import Base: endof, length

using IndexableDicts
import IndexableDicts: rank

immutable WaveletMatrix{B<:AbstractBitVector}
    bits::Vector{B}
    nzeros::Vector{Int}
    len::Int
end

function WaveletMatrix(data::Vector{Uint8})
    n = length(data)
    bits = SuccinctBitVector[]
    nzeros = Int[]
    # TODO: efficient construction
    data′ = Array(Uint8, n)
    for d in 1:8
        # scan d-th bit
        bits′ = SuccinctBitVector()
        for i in 1:n
            if bitat(data[i], 8 - d + 1)
                # right
                push!(bits′, 1)
            else
                # left
                push!(bits′, 0)
            end
        end
        nzero = rank0(bits′, n)
        l = r = 1
        for i in 1:n
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
    WaveletMatrix(bits, nzeros, n)
end

WaveletMatrix(str::ByteString) = WaveletMatrix(copy(str.data))

endof(wm::WaveletMatrix) = wm.len
length(wm::WaveletMatrix) = wm.len

function rank(a::Uint8, wm::WaveletMatrix, i::Int)
    sp = 0
    ep = i
    # scan from the most significant bit to the least significant bit
    # e.g.
    #   0b01000101
    #     =======>
    #   d=1 .... 8
    for d in 1:8
        bit = bitat(a, 8 - d + 1)
        #@show bit sp:ep
        sp = rank(bit, wm.bits[d], sp)
        ep = rank(bit, wm.bits[d], ep)
        if bit
            sp += wm.nzeros[d]
            ep += wm.nzeros[d]
        end
        #@show sp:ep
    end
    return ep - sp
end

rank(c::Char, wm::WaveletMatrix, i::Int) = rank(convert(Uint8, c), wm, i)

function bitat(byte::Uint8, i::Int)
    return ((byte >> (i-1)) & 0x01) == 1
end

end # module
