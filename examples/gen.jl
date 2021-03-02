include("../src/MatrixGen.jl")
using MatrixMarket

folder = "mats"
@show folder 

d = 3
for n = [10,20,30,40,50,60]
    for rho = [1, 10, 100, 1000]
        @show n, rho
        (A, a) = Ad_hc(n, d, rho)
        mmwrite("$(folder)/neglapl_$(d)_$(n)_$(rho).mm", A)
    end
end
