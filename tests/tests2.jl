include("../src/MatrixGen.jl")
A, X = elliptic_dirichlet(3, 2, x -> 1, x -> [1, 0], x -> 0)
@show A
@show X