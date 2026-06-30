#!/usr/bin/env julia
# Convert the Quaia-constrained white noise to a MUSIC white-noise file, read per level via
#     [random]  seed[<level>] = <this file>
#
#   julia make_music_wnoise.jl <carrier.npz> <res> <out.bin> [seed]
#
# where 2^level == res (e.g. res=1024 → level 10, res=256 → level 8). Point MUSIC's seed[level] at the
# produced file (see examples/music_unigrid.conf and music_zoom.conf).
#
# MUSIC white-noise format (Fortran-unformatted, 32-bit record markers), reverse-engineered from
# MUSIC's plugins/random_music_wnoise_generator.cc:
#   header record :  [int32 16] [uint32 nx] [uint32 ny] [uint32 nz] [int32 iseed] [int32 16]
#   then nz records (one per z-plane k=1..nz), each:
#                    [int32 nx*ny*4] [float32 × nx*ny, x fastest then y] [int32 nx*ny*4]
# The (x fastest, then y, planes stacked in z) ordering is Julia's native column-major layout, so each
# z-slice ω[:, :, k] writes directly. MUSIC negates the field on read by default; set
#     [random]  grafic_sign = yes
# so it uses the field as-is. VERIFY the sign on first use — overdensities should land where the Quaia
# quasars are (compare a smoothed slice to the catalog).

using DiscoInverse, NPZ

carrier = ARGS[1]
res     = parse(Int, ARGS[2])
out     = ARGS[3]
seed    = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 101
ispow2(res) || @warn "res=$res is not a power of two; MUSIC levels are 2^n"

d  = npzread(carrier)
ω  = Float32.(refine_phases(d["phi_coarse"], res; seed=seed))
nx, ny, nz = size(ω)
open(out, "w") do io
    write(io, Int32(16), UInt32(nx), UInt32(ny), UInt32(nz), Int32(seed), Int32(16))   # header record
    bs = Int32(nx * ny * sizeof(Float32))
    for k in 1:nz
        write(io, bs); write(io, @view ω[:, :, k]); write(io, bs)                       # z-plane record
    end
end
println("wrote MUSIC white noise ", out, "  (", nx, "^3 → level ", round(Int, log2(res)),
        ")  — set seed[", round(Int, log2(res)), "] = ", out, " and grafic_sign = yes")
