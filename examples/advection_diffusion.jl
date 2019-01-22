using Plots

## Example enerating advection-diffusion matrices and solution
include("../src/MatrixGen.jl")

n = 64
d = 2
A, X = elliptic_dirichlet(n, d, x -> 1, x -> [30, 30], x -> 0)
f = rand(size(A,1))
u = A\f

U = reshape(u, (n, n))
F = reshape(f, (n, n))

p1 = contourf(U)
display(p1)