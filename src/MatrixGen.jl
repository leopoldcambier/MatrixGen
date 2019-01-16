using ImageFiltering
using SparseArrays
using Random
using LinearAlgebra

## From map to sparse CSC matrix
# Converts a map (i,j)=>v to a sparse matrix, where i, j are D-tuples of integers
# In:
#   - A: a map (i,j) => v
# Out:
#   - I: I[k] = i indicates tuple i corresponds to row/column k in A
#   - A: A[i,j] = v, the sparse matrix
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

# Creates a linear system corresponding to 
# - div(a(x) grad(u(x))) + b(x) grad(u(x)) + c(x) u(x) = f(x)     (*)
# with Dirichlet boundary conditions over [0,1]^d, 
# discretized using n+2 points in each dimension. 
# The boundaries are set to 0 and, hence, not included
# b(x) is assumed to be incompressible (div(b(x)) = 0 for all x).
#
# In each dimension, we have
# U(i-1) * ( -a(i-1/2)/h^2 - b(i-1)[d]/2h      )
# U(i)   * ( (a(i-1/2)+a(i+1/2))/h^2 + c(i) )
# U(i+1) * ( -a(i+1/2)/h^2 + b(i+1)[d]/2h      )
#
# The resulting matrix is of size n^d x n^d
# Returns (A, X) where X is d x n^d, the coordinates of each point on the grid and
# A is n^d x n^d sparse matrix#
#
# In:
#   - n: number of discretization points in each dimension
#   - d: number of dimensions
#   - a: function [0,1]^d -> R, the diffusion coefficient
#   - b: function [0,1]^d -> R^d, the flow field
#   - c: function [0,1]^d -> R
# Out:
#   - A: the n^d x n^d sparse CSC matrix corresponding to the discretization of (*) with second-order FD
#   - X: the d x n^d grid coordinates of every point
function elliptic_dirichlet(n, d, a, b, c)
    A = Dict{Tuple{NTuple{d,Int},NTuple{d,Int}},Float64}()
    h = 1/(n+1)
    for i_ in CartesianIndices(tuple([1:n for i in 1:d]...))
        i = Tuple(i_)
        for dim = 1:d
            x = h .* i
            e = Tuple([Int(i == dim) for i = 1:d])
            A[(i,i .- e)] =                  -  a(x .- h.*e./2)/h^2                   - b(x .- h.*e)[dim]/(2h)
            A[(i,i     )] = get(A, (i,i), 0) + (a(x .- h.*e./2) + a(x .+ h.*e./2))/h^2                         + c(x)
            A[(i,i .+ e)] =                  -  a(x .+ h.*e./2)/h^2                   + b(x .+ h.*e)[dim]/(2h)
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

# Estimate the condition number of an SPD matrix A
# In: A, a sparse SPD matrix
# Out: An estimate for cond(A)
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
# In: 
#   - n: the number of vertices
#   - deg: the degree of each vertex
# Out:
#   - A: a sparse adjacency matrix of the graph
#   - x: coordinate of vertices
#   - y: coordinate of vertices
#   - z: coordinate of vertices
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

# Generate high-contrast semi-random images
# In:
#   - n: each dimension size
#   - d: number of dimension
#   - rho: highs are rho and lows are 1/rho
#   - sigma: typical size of the patches
# Out:
#   - a: n x n ... x n tensor
function high_contrast_field(n, d, rho, sigma)
    a = rand(MersenneTwister(0), [n for i = 1:d]...)    # Uniformly random n x n ... x n trensor
    Ig = imfilter(a, Kernel.gaussian(sigma))            # Smooth
    If = (Ig .<= 0.5)                                   # Quantize
    a[If] .= rho                                        # Highs
    a[.! If] .= 1/rho                                   # Lows
    return a
end

# A high-contrast laplacian
function Ad_hc(n, d, rho)

    # Build quantized high contrast field
    a = high_contrast_field(2n+1, d, rho, 2)

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

