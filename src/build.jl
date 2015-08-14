using IntArrays

# build an internal structure of the wavelet matrix
function build{T<:Unsigned,B}(::Type{B}, data::AbstractVector{T}, n::Int)
    ivec = IntVector{n}(data)
    bits = B[]
    _build!(B, ivec, bits, n)
    return tuple(bits...)
end

function _build!{B}(::Type{B}, ivec, bits, n)
    len = length(ivec)
    # allocate working space
    bv = BitVector(len)
    ivec′ = similar(ivec)
    for d in 1:n
        # scan d-th bit
        nzeros = 0
        for i in 1:len
            if (ivec[i] >> (n - d)) & 1 == 1
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
        push!(bits, B(bv))
        copy!(ivec, ivec′)
    end
end

# starting points for the rank operation can be precomputed
function locate_sp(a, bits, nzeros)
    n = length(bits)
    sp = 0
    for d in 1:n
        bit = (a >> (n - d)) & 1 == 1
        sp = rank(bit, bits[d], sp)
        if bit
            sp += nzeros[d]
        end
    end
    return sp
end
