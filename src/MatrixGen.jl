using ImageFiltering
using SparseArrays
using Random

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
    Fc = cholfact(A)
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
