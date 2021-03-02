include("../src/MatrixGen.jl")
using Plots
gr()

## Generate high contrast laplacians
for n = [50,63,80,101,127,160,202,254]
    for rho = [1,100]
        println("Rho = ", rho, " n = ", n)
        (A, a) = Ad_hc(n, 3, rho)
        println("Cond(A) +-= ", estimate_cond2(A))
        heatmap(a, aspect_ratio=1)
    end
end
#heatmap(a, aspect_ratio=1)