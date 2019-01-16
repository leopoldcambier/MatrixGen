using Plots
pyplot()

## Example enerating advection-diffusion matrices and solution
include("../src/MatrixGen.jl")

# Generate an image with patches
n = 128
d = 2
rho = 10
sigma = 10
a = high_contrast_field(n+2, d, rho, sigma)
# Define a(x) 
function af(x)
    i = Int.(round.((n+1) .* x)).+1
    return a[i...]
end

# Generate matrix
A, _ = elliptic_dirichlet(n, d, af, x -> [30, 30], x -> 0)

# Random rhs
f = rand(size(A,1))

# Solve
u = A\f

# Plot
U = reshape(u, (n, n))
F = reshape(f, (n, n))
display(contourf(U, reuse=false))
display(contourf(F, reuse=false))
display(contourf(a, reuse=false))