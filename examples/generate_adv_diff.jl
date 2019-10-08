include("../src/MatrixGen.jl")
using MatrixMarket
using Printf
using ArgParse

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin
        "--n"
            help = "n, the 1D size"
            arg_type = Int
        "--d"
            help = "d, the number of dimensions"
            arg_type = Int
        "--a"
            help = "a, the diffusion coefficient"
            arg_type = Float64
        "--b"
            help = "b, the advection coefficient"
            arg_type = Float64
        "--dir"
            help = "the file prefix (directory)"
            arg_type = String
    end
    return parse_args(s)
end

function generate_advdiff_ab(n, d, a, b, directory)
    file = string(directory, @sprintf("/advdiff_centered_const_%d_%d_%d_%d.mm", d, n, a, b))
    @printf("Generating advdiff problem with d = %d, n = %d, a = %f, b = %f, in %s\n", d, n, a, b, file)
    @time A, _ = elliptic_dirichlet(n, d, x -> a, x -> [b for i= 1:d], x -> 0, upwind=false)
    @time mmwrite(file, A)
    @printf("Done generating")
end

function main()
    pa = parse_commandline()
    println(pa)
    generate_advdiff_ab(pa["n"], pa["d"], pa["a"], pa["b"], pa["dir"])
end

main()
