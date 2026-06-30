#!/usr/bin/env julia
# Realize the Quaia-constrained white-noise box at ANY resolution from the constraint carrier.
#
#   julia realize_box.jl <carrier.npz> <res_box> <out.f32> [seed]
#
# The carrier holds the constrained coarse phases (res = 256). `refine_phases` embeds them into the
# fine grid at their physical wavenumbers and fills the new small-scale modes with fresh random phases
# — a pure FFT (CPU, no GPU). The result is a unit-variance periodic white-noise field with exactly the
# fiducial P(k); its large scales carry the Quaia constraint, its small scales are a fresh ΛCDM draw.
#
# Requires DiscoInverse.jl (which pulls in DiscoDJ.jl):
#   julia> import Pkg; Pkg.add(url="https://github.com/yipihey/DiscoInverse.jl")

using DiscoInverse, NPZ

carrier = ARGS[1]
res_box = parse(Int, ARGS[2])
out     = ARGS[3]
seed    = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 101

d  = npzread(carrier)
println("carrier: res_constrain=", Int(d["res_constrain"]),
        "  L_box=", round(d["boxsize"]), " Mpc/h  constrained_radius=", round(d["constrained_radius"]), " Mpc/h")

ω = refine_phases(d["phi_coarse"], res_box; seed=seed)   # (res_box, res_box, res_box), unit-variance white
open(io -> write(io, Float32.(ω)), out, "w")             # raw little-endian float32, column-major (x fastest)
println("wrote ", out, "  (", res_box, "^3 float32, ",
        round(res_box^3 * 4 / 1e9; digits=2), " GB, seed=", seed, ")")
