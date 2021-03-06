# ACRM.jl

This code base is using the Julia Language and [DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> ACRM.jl

It is authored by Reza Arefidamghani, Roger Behling, Alfredo N. Iusem, Luiz-Rafael Santos.

To (locally) reproduce this project, do the following:

0. Download this code base. Notice that raw data are typically not included in the
   git-history and may need to be downloaded independently. The preferred way would be
   to clone the repo using:
   ```
   git clone https://github.com/lrsantos11/ACRM.jl.git
   ```
   or better, fork to you GitHub page and clone from there.
1. Open a Julia console and do:
   ```julia
   julia> using Pkg
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box.
