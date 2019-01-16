using ImageFiltering
using SparseArrays
using Random
using LinearAlgebra

# Return e_k vector
function ed(d,k)
    x = zeros(Int64,d)
    x[k] = 1
    return tuple(x...)
end

# Convert a map{(i,j)=>v} into a sparse matrix
# The i, j's have to be sortable
function map2mat(A::Dict{Tuple{NTuple{D,Int},NTuple{D,Int}},Float64}) where D
    I = sort(unique([i for (i,j) in keys(A)]))
    J = sort(unique([j for (i,j) in keys(A)]))
    @assert(all(I .== J))
    map = Dict{NTuple{D,Int},Int}()
    for i = 1:length(I)
        map[I[i]] = i
    end
    II = I
    m = length(II)
    n = length(A)
    I = zeros(Int64, n)
    J = zeros(Int64, n)
    V = zeros(Float64, n)
    for (k,((i,j),v)) in enumerate(A)
        I[k] = map[i]
        J[k] = map[j]
        V[k] = v
    end
    return II, sparse(I,J,V,m,m) # Ordering, Matrix
end

# In each dimension, we have
# U(i-1) * ( -a(i-1/2)/h^2 - b(i-1)[d]/2h      )
# U(i)   * ( (a(i-1/2)+a(i+1/2))/h^2 + c(i) )
# U(i+1) * ( -a(i+1/2)/h^2 + b(i+1)[d]/2h      )
# We have n+2 points over [0, 1] but boundaries = 0 and are excluded
# The resulting matrix is of size n^d x n^d
# Returns (A, X) where X is d x n^d, the coordinates of each point on the grid
# a : R^d -> R, a(x) > 0
# b : R^d -> R^d
# c : R^d -> R
function elliptic_dirichlet(n, d, a, b, c)
    A = Dict{Tuple{NTuple{d,Int},NTuple{d,Int}},Float64}()
    h = 1/(n+1)
    for i_ in CartesianIndices(tuple([1:n for i in 1:d]...))
        i = Tuple(i_)
        for dim = 1:d
            x = h .* i
            e = Tuple([Int(i == dim) for i = 1:d])
            A[(i,i .- e)] =                  -  a(x .- e./2)/h^2                 - b(x .- e)[dim]/(2h)
            A[(i,i     )] = get(A, (i,i), 0) + (a(x .- e./2) + a(x .+ e./2))/h^2                       + c(x)
            A[(i,i .+ e)] =                  -  a(x .+ e./2)/h^2                 + b(x .+ e)[dim]/(2h)
        end
    end
    inside = x -> all((x .>= 1) .& (x .<= n))
    A = filter(x -> inside(x.first[1]) && inside(x.first[2]), A)
    I, A = map2mat(A)
    X = zeros(Float64, (d, n^d))
    @assert(length(I) == n^d)
    for (k,i) = enumerate(I)
        X[:,k] = h * [i...]
    end
    return A, X
end

# Generate a matrix [dim x N]
# where each column are some coordinates
# in a tensor [x1 * x2 * ... * xd] space
function tensor_grid(grids)
    n = 1
    for g in grids
        n *= length(g)
    end
    X = zeros(length(grids), n)
    for (i,x) in enumerate(Base.product(grids...))
        X[:,i] = [x...]
    end
    return X
end

# A [-1, 2, -1] 3-points stencil matrix
# with Dirichlet BC's
An = n -> spdiagm(-1 => -ones(n-1), 0 => 2*ones(n),1 => -ones(n-1))

# A 3, 5, or 7-points stencil negative Laplacian
# on n^d cube with Dirichlet BC's
function Ad(n, d)
    I = spdiagm(0 => ones(n))
    if d == 1
        A = An(n)
    elseif d == 2
        A = kron(An(n), I) + kron(I, An(n))
    elseif d == 3
        A = kron(An(n), I, I) + kron(I, An(n), I) + kron(I, I, An(n))
    end
    return A
end

# A high-contrast laplacian
function Ad_hc(n, d, rho)

    # Build quantized high contrast field
    @assert(rho >= 1)
    a = rand(MersenneTwister(0), [2n+1 for i = 1:d]...)
    Ig = imfilter(a, Kernel.gaussian(2)) # 2 since 2*n+1
    If = (Ig .<= 0.5)
    a[If] .= rho
    a[.! If] .= 1/rho

    # Assemble matrix
    I = zeros(Int64, 0)
    J = zeros(Int64, 0)
    V = zeros(Float64, 0)

    sizehint!(I, (2*d+1)*n)
    sizehint!(J, (2*d+1)*n)
    sizehint!(V, (2*d+1)*n)
    
    @assert d <= 3

    if d == 3
        L  = [ (-1, 0, 0), (1, 0, 0), (0, -1, 0), (0, 1, 0), (0, 0, -1), (0, 0, 1) ]
        idxf = (i,j,k) -> i + (j - 1) .* n + (k - 1) .* n .* n 
        for i = 1:n
            for j = 1:n
                for k = 1:n
                    idx = idxf(i,j,k)
                    av = [a[2*i+i_,
                            2*j+j_,
                            2*k+k_] for (i_,j_,k_) in L]
                    push!(I, idx)
                    push!(J, idx)
                    push!(V, sum(av))
                    for ((i_,j_,k_), av_) in zip(L, av)
                        idx_ = idxf(i+i_,j+j_,k+k_)
                        if (i+i_ >= 1 && i+i_ <= n &&
                            j+j_ >= 1 && j+j_ <= n &&
                            k+k_ >= 1 && k+k_ <= n)
                            push!(I, idx)
                            push!(J, idx_)
                            push!(V, - av_)
                        end
                    end
                end
            end
        end
    elseif d == 2
        L =  [ (-1, 0), (1, 0), (0, -1), (0, 1) ]
        idxf = (i,j) -> i + (j - 1) .* n
        for j = 1:n
            for i = 1:n
                idx = idxf(i,j)
                av = [a[2*i+i_,
                        2*j+j_] for (i_,j_) in L]
                push!(I, idx)
                push!(J, idx)
                push!(V, sum(av))
                for ((i_,j_), av_) in zip(L, av)
                    idx_ = idxf(i+i_,j+j_)
                    if (i+i_ >= 1 && i+i_ <= n &&
                        j+j_ >= 1 && j+j_ <= n)
                        push!(I, idx)
                        push!(J, idx_)
                        push!(V, - av_)
                    end
                end
            end
        end
    end

    A = sparse(I, J, V, n^d, n^d)

    return (A, a)
end

# Good but slow
function estimate_cond2(A)

    nA, nAi = 0, 0
    Fc = cholesky(Hermitian(A))
    doneA = false
    doneAi = false

    x = randn(size(A,1))
    n = 200

    Anx = copy(x)
    for i = 1:n
        Anx_ = A*Anx
        nA_ = norm(Anx_) / norm(Anx)
        if abs(nA - nA_) / nA_ < 1e-3
            doneA = true
            println("||A|| converged in $i iterations")
            break
        end
        nA = max(nA_, nA)
        Anx = Anx_
    end

    Ainx = copy(x)
    for i = 1:n
        Ainx_ = Fc\Ainx
        nAi_ = norm(Ainx_) / norm(Ainx)
        if abs(nAi - nAi_) / nAi_ < 1e-3
            doneAi = true
            println("||A^-1|| converged in $i iterations")
            break
        end
        nAi = max(nAi, nAi_)
        Ainx = Ainx_
    end

    if (!doneA) || (!doneAi)
        println("Reached max number of iterations in cond2 estimate")
    end

    return nA * nAi

end

# Random graph in 3D space with connections
# to nearest neighbors
function rand_A_3D(n::Int64, deg::Int64) 
    @assert deg < n
    XYZ = rand(n, 3)
    x = XYZ[:,1]
    y = XYZ[:,2]
    z = XYZ[:,3]
    I, J, V = Int64[], Int64[], Float64[]
    sizehint!(I, deg*n)
    sizehint!(J, deg*n)
    sizehint!(V, deg*n)
    for i = 1:n
        d = (x-x[i]).^2 + (y-y[i]).^2 + (z-z[i]).^2
        p = sortperm(d)
        jj = p[2:(2+deg-1)]
        jj = jj[jj .> i]
        push!(I, i)
        push!(J, i)
        push!(V, deg)
        for j in jj
            push!(I, i)
            push!(J, j)
            push!(V, -1)
        end
    end
    A = sparse(I, J, V, n, n)
    A = (A + A')/2.0
    return (A, x, y, z)
end
