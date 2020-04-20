include("../src/MatrixGen.jl")
using Plots

## Generate high contrast laplacians
for n = [200, 400, 800, 1600, 3200, 6400]
    for rho = [100, 1000]
        println("Rho = ", rho, " n = ", n)
        (A, a) = Ad_hc(n, 2, rho)
        println("Cond(A) +-= ", estimate_cond2(A))
    end
end
#heatmap(a, aspect_ratio=1)
