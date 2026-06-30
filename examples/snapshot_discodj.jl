#!/usr/bin/env julia
# Make a particle IC snapshot directly from the carrier with DiscoDJ.jl — no external IC code.
#
#   julia snapshot_discodj.jl <carrier.npz> <res> <z_init> [seed]
#
# Realizes the white noise at `res`, applies the IC operator + 2LPT at scale factor a = 1/(1+z_init),
# and returns periodic-wrapped Lagrangian-displaced particle positions. This is the same forward DiscoDJ
# uses internally, so the constrained large-scale field is preserved exactly. Saves positions to .npz.
#
# Note: this writes positions only (a demonstration). For a ready-to-run Gadget/Enzo file, either route
# the white noise through MUSIC (examples/make_music_wnoise.jl + music_unigrid.conf), or add 2LPT
# velocities (v ∝ a²Ḋ ψ) and your sim's header with DiscoDJ's writers.

using DiscoInverse, NPZ

carrier = ARGS[1]
res     = parse(Int, ARGS[2])
zinit   = parse(Float64, ARGS[3])
seed    = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 101

d   = npzread(carrier); L = d["boxsize"]; c = fiducial_cosmology()
ω   = refine_phases(d["phi_coarse"], res; seed=seed)
pos = ic_box_snapshot(ω, L, c, 1 / (1 + zinit); n_order=2)        # 2LPT positions, (res,res,res,3), periodic
outp = replace(carrier, ".npz" => string("_snap_z", zinit, "_", res, ".npz"))
npzwrite(outp, Dict("pos" => pos, "boxsize" => L, "z_init" => zinit, "res" => res))
println("2LPT snapshot at z=", zinit, ": ", size(pos), " positions in [0,", round(L), ") Mpc/h → ", outp)
