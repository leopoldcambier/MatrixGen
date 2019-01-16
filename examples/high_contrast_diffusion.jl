include("../src/MatrixGen.jl")
using Plots

## Generate high contrast laplacians
for n = [512, 1024, 2048]
    for rho = [1, 10, 100, 1000]
        println("Rho = ", rho, " n = ", n)
        (A, a) = Ad_hc(n, 2, rho)
        println("Cond(A) +-= ", estimate_cond2(A))
    end
end
#heatmap(a, aspect_ratio=1)