using IntArrays

# build an internal structure of the wavelet matrix
function build(::Type{B}, data::AbstractVector{T}, w::Int, destructive::Bool) where {T<:Unsigned,B}
    ivec = IntVector{w}(data)
    if destructive
        # free memory
        empty!(data)
    end
    bits = Vector{B}(undef,w)
    _build!(B, ivec, bits)
    return tuple(bits...)
end

function _build!(::Type{B}, ivec, bits) where B
    w = length(bits)
    len = length(ivec)
    # allocate working space
    bv = BitVector(undef, len)
    ivec′ = similar(ivec)
    for d in 1:w
        # scan d-th bit
        nzeros = 0
        for i in 1:len
            if (ivec[i] >> (w - d)) & 1 == 1
                bv[i] = true  # right
            else
                bv[i] = false # left
                nzeros += 1
            end
        end
        # stably sort integers by bv
        r = nzeros
        l = 0
        for i in 1:len
            if bv[i]
                ivec′[r+=1] = ivec[i]
            else
                ivec′[l+=1] = ivec[i]
            end
        end
        # store the bit vector and go next
        bits[d] = bv
        copy!(ivec, ivec′)
    end
end

# starting points for the rank operation can be precomputed
function locate_sp(a, bits, nzeros)
    w = length(bits)
    sp = 0
    for d in 1:w
        bit = (a >> (w - d)) & 1 == 1
        sp = rank(bit, bits[d], sp)
        if bit
            sp += nzeros[d]
        end
    end
    return sp
end
