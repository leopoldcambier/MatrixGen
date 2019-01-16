include("../src/MatrixGen.jl")
using MatrixMarket
using Printf

for (d, n) in [ (2, 3), (2, 5), (2, 10), (2, 16), (2, 20), (2, 32), (2, 64), (2, 96), (2, 100), (2, 128), (2, 256), 
    (3, 5), (3, 10), (3, 15), (3, 16), (3, 25), (3, 30), (3, 32), (3, 48), (3, 64), (3, 96)]
    @show d, n
    A = Ad(n, d)
    (I,J,V) = findnz(A)
    V += 0.1 * randn(size(V))
    A = sparse(I,J,V)
    mmwrite(@sprintf("neglapl_unsym_%d_%d.mm", d, n), A)
end

