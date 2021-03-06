include("../src/MatrixGen.jl")
using Plots
gr()

## Generate high contrast laplacians
for n = [30]
    for rho = [100]
        println("Rho = ", rho, " n = ", n)
        (A, a) = Ad_hc(n, 2, rho)
        println("Cond(A) +-= ", estimate_cond2(A))
        heatmap(a, aspect_ratio=1)
    end
end
#heatmap(a, aspect_ratio=1)