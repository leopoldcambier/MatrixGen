include("../src/MatrixGen.jl")
using MatrixMarket
using Printf

for b = [10, 10000]
    for (d, n) in [ (2, 256), (2, 512), (2, 1024), (2, 2048), (2, 4096)]
        @show d, n
        h = 1/(n+1)
        A, _ = elliptic_dirichlet(n, d, x -> 1, x -> [b for i= 1:d], x -> 0)
        mmwrite(@sprintf("advdiff_%d_%d_1_%d.mm", d, n, b), A)
    end
end
