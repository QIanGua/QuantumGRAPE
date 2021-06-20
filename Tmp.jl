module Tmp

include("vars.jl")
# include("main.jl")

# export functions
export say_hello
say_hello() = println("Hello!")



function Hilbert_to_Liouville_index(n::Int64)
    index = zeros(Int64,n,n)
    step = 1
    for i = 1:n
        index[ i,i ] = step
        step += 1
    end

    for i = 1:n
        for j = i+1:n
            index[ i,j ] = step
            step += 1
            index[ j,i ] = step
            step += 1
        end
    end
    return index
end

function CommutatorMatrix(H)
    n = size(H,1)
    L = zeros(n^2,n^2)
    index = Hilbert_to_Liouville_index(n)
    for m1 = 1:n^2
        for n1 = 1:n^2
            i1 = findall(x->x==m1,index)[1][1]
            j1 = findall(x->x==m1,index)[1][2]
            i2 = findall(x->x==n1,index)[1][1]
            j2 = findall(x->x==n1,index)[1][2]

            L[m1,n1] = H[i1,i2]*(j1==j2) - (i1 == i2)*H[j2,j1]
        end
    end
    return L
end

function AntiCommutatorMatrix(H)
    n = size(H,1)
    L = zeros(n^2,n^2)
    index = Hilbert_to_Liouville_index(n)
    for m1 = 1:n^2
        for n1 = 1:n^2
            i1 = findall(x->x==m1,index)[1][1]
            j1 = findall(x->x==m1,index)[1][2]
            i2 = findall(x->x==n1,index)[1][1]
            j2 = findall(x->x==n1,index)[1][2]
            L[m1,n1] = H[i1,i2]*(j1==j2) + (i1 == i2)*H[j2,j1]
        end
    end
    return L
end

end