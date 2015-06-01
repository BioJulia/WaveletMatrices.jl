# FM-index searching for ASCIIString
# Reference: PAOLO, et al. "Compressed Text Indexes: From Theory to Practice"

# https://github.com/quinnj/SuffixArrays.jl
using SuffixArrays
using WaveletMatrices
import WaveletMatrices: rank

# Return the interval SA[sp,ep] of text suffixes prefixed by p; see "Algorithm FM-count"
# Arguments:
#   t: BWTed text
#   p: prefix
#   c: character counts
function countFM(t, p, cnt)
    sp, ep = 1, endof(t)
    i = endof(p)
    while sp ≤ ep && i ≥ 1
        c = p[i]
        sp = cnt[c] + rank(c, t, sp - 1) + 1
        ep = cnt[c] + rank(c, t, ep)
        i -= 1
    end
    return ifelse(sp > ep, 0:-1, sp:ep)
end

function bwt(t)
    sa = suffixsort(t)
    # note that the indices stored in a saffix array is 0 based
    ord = [sa.index[i] > 0 ? sa.index[i] : endof(t) for i in 1:endof(t)]
    return t[ord], sa
end

function count(t)
    cnt = zeros(Int, 128)
    for c in t
        cnt[convert(Int, c)+1] += 1
    end
    return cumsum(cnt)
end

# poor man's rank query (should be implemented in a Wavelet tree/matrix)
function rank(c, t, i)
    @assert 0 <= i <= endof(t)
    n = 0
    for i′ in 1:i
        if t[i′] == c
            n += 1
        end
    end
    return n
end

function count_search(text, query)
    n = 0
    start = 1
    while true
        loc = search(text, query, start)
        if isempty(loc)
            return n
        end
        start = first(loc) + 1
        n += 1
    end
end

function main()
    textfile = ARGS[1]

    # build a self-indexed text
    info("loading a text")
    text = replace(readall(textfile), '\n', "")
    info("building a BWT")
    t, _ = bwt(text)
    info("counting characters")
    cnt = count(t)
    info("building a WaveletMatrix")
    wm = WaveletMatrix(t)

    # read, print, search, loop
    print("query> ")
    query = readline(STDIN)
    while !isempty(query)
        query = chomp(query)
        if query[end] == ';'
            query = chop(query)
            gc()
            t0 = time_ns()
            n = count_search(text, query)
        else
            gc()
            t0 = time_ns()
            n = length(countFM(wm, query, cnt))
        end
        t1 = time_ns()
        println("'$query' occurs $n times ($(div(t1 - t0, 1000)) μs).")
        print("query> ")
        query = readline(STDIN)
    end
end

main()
