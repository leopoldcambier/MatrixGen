# MatrixGen

Basic functions to generate basic problems (laplacians, high-contrast laplacians, simple advection-diffusion, etc)
to benchmark linear solvers

Still a very very basic code

For instance, you can generate the corresponding images
![a](https://github.com/leopoldcambier/MatrixGen/raw/master/pics/a.png)

Then you solve an advection-diffusion using it as a background field, leading for instance to
![u](https://github.com/leopoldcambier/MatrixGen/raw/master/pics/u.png)

See examples/advection_diffusion_high_contrast.jl for an example