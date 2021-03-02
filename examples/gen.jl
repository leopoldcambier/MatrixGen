include("../src/MatrixGen.jl")
using MatrixMarket

folder = "mats"
@show folder 

d = 3
for n = [50,63,80,101,127,160,202,254]
    for rho = [1,100]
        @show n, rho
        (A, a) = Ad_hc(n, d, rho)
        mmwrite("$(folder)/neglapl_$(d)_$(n)_$(rho).mm", A)
    end
end
