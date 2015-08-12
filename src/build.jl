using IntArrays

# build an internal structure of the wavelet matrix
function build{T<:Unsigned,B}(::Type{B}, data::AbstractVector{T}, n::Int)
    len = length(data)
    bits = B[]
    ivec = IntVector{n}(data)
    for d in 1:n
        # scan d-th bit
        bv = BitVector(len)
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
        ivec′ = similar(ivec)
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
        ivec = ivec′
    end
    return tuple(bits...)
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
