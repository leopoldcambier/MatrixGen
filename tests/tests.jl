include("../src/MatrixGen.jl")
using Plots

# 2-D Laplacian
A = Ad(5, 2)

# High contrast
(A, a) = Ad_hc(100, 2, 2)

heatmap(a)
