using Plots

include("../src/MatrixGen.jl")
n = 64
d = 2
A, X = elliptic_dirichlet(n, d, x -> 1, x -> [0, 0], x -> 0)
f = [sum(X[:,i]) for i = 1:size(X,2)]
u = A\f
U = reshape(u, (n, n))
p2 = contourf(U)