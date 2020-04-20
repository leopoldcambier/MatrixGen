include("../src/MatrixGen.jl")
using MatrixMarket

folder = "mats"
@show folder 

d = 2
for n = [100,200,400,800,1600,3200,6400]
    for rho = [1, 10, 100, 1000, 10000]
        @show n, rho
        (A, a) = Ad_hc(n, d, rho)
        mmwrite("$(folder)/neglapl_$(d)_$(n)_$(rho).mm", A)
    end
end
