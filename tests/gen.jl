include("../src/MatrixGen.jl")
using MatrixMarket

folder = "/scratch/groups/darve/lcambier/matrices/laplacians2/"

d = 2
for n = [512, 1024, 2048, 4096, 8192]
    for rho = [1, 10, 100, 1000, 10000]
        (A, a) = Ad_hc(n, d, rho)
        mmwrite("$(folder)/neglapl_$(d)_$(n)_$(rho).mm", A)
    end
end
