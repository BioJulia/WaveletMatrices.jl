# build an internal structure of the wavelet matrix
function build{T,B}(::Type{B}, data::Vector{T}, n::Int)
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
    return tuple(bits...), tuple(nzeros...)
end
