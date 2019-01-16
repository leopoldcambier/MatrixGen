include("../src/MatrixGen.jl")
using MatrixMarket
using Printf

for (d, n) in [ (2, 256), (2, 512), (2, 1024), (3, 16), (3, 32), (3, 64), (3, 96), (3, 128) ]
    @show d, n
    h = 1/(n+1)
    a = high_contrast_field(n+2, d, 10, 5)
    function af(x)
        i = Int.(round.((n+1) .* x)).+1
        return a[i...]
    end
    A, _ = elliptic_dirichlet(n, d, af, x -> [10 for i= 1:d], x -> 0)
    mmwrite(@sprintf("advdiff_%d_%d.mm", d, n), A)
end

