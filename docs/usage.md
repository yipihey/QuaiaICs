# Usage — running simulations from the Quaia-constrained ICs

## 0. Get the carrier and the tools

- Download the constraint carrier `quaia_icbox_phases.npz` (~67 MB) — see [../data/README.md](../data/README.md).
- Install [DiscoInverse.jl](https://github.com/yipihey/DiscoInverse.jl) (pulls in DiscoDJ.jl):
  ```julia
  import Pkg; Pkg.add(url="https://github.com/yipihey/DiscoInverse.jl")
  ```
  Only `refine_phases`, `fiducial_cosmology`, `ic_box_snapshot` and `NPZ` are needed to *use* the ICs —
  no GPU, the realization is a CPU FFT.

## 1. Realize the white-noise box at your resolution

```bash
julia examples/realize_box.jl quaia_icbox_phases.npz 1024 omega_1024.f32        # 1024^3
julia examples/realize_box.jl quaia_icbox_phases.npz 2048 omega_2048.f32 7      # 2048^3, seed 7
```

`omega_*.f32` is a raw little-endian `float32`, column-major (x fastest), `res^3` contiguous — a
unit-variance periodic white-noise field. The large-scale phases inside the Quaia sphere are constrained;
the small scales and the box corners are a fresh ΛCDM draw.

## 2. Choosing resolution and box size

- **Resolution = your choice.** The constraint lives at k ≲ 0.05 h/Mpc; any grid whose Nyquist exceeds
  that resolves it. The small scales are filled by `refine_phases` with a statistically-correct ΛCDM
  random field that inherits the constrained large-scale tidal field. Common choices: 512³ (testing),
  1024³–2048³ (production).
- **Box size** for *this* carrier is fixed at L = 13 260 Mpc/h (it must contain the all-sky Quaia
  volume). To realize into a **larger** periodic box (more unconstrained padding around the Quaia
  sphere), regenerate a carrier with a bigger `boxsize` in DiscoInverse:
  ```julia
  using DiscoInverse
  box = constrained_ic_box(cat, randoms, cosmo; L_box=20000, res_constrain=256, res_box=1024)
  export_white_noise("box.f32", box)            # field + box.f32.manifest.npz
  ```
  A *smaller* box would cut the survey and is not supported.
- **Ensemble:** different `seed`s give independent realizations sharing the constraint — generate as many
  as you need for an uncertainty ensemble.

## 3. MUSIC (recommended bridge to any sim code)

MUSIC reads an external white-noise field per level and outputs Gadget/Enzo/grafic ICs, with full
uni-grid and zoom support.

**Uni-grid** (`examples/music_unigrid.conf`):
```bash
julia examples/make_music_wnoise.jl quaia_icbox_phases.npz 1024 wnoise_music_1024.bin
MUSIC examples/music_unigrid.conf        # → ic_quaia_1024.gdt
```
Key points: `2^levelmax == resolution` of the white-noise file (1024 → level 10); `boxlength` must equal
the carrier `boxsize` (13260); set `grafic_sign = yes` so MUSIC uses the field as-is; pick the output
`format` (gadget2 / enzo / grafic / generic). **Verify the sign on first use** — a smoothed slice of the
resulting δ should have overdensities where the Quaia quasars cluster.

**Zoom-in** (`examples/music_zoom.conf`): give MUSIC the *constrained coarse field* at the base level and
let MUSIC's own multi-scale white-noise refinement (Hahn & Abel 2011) add small-scale power inside the
high-resolution Lagrangian patch. The Quaia-constrained modes then provide the correct large-scale tidal
field around the zoom:
```bash
julia examples/make_music_wnoise.jl quaia_icbox_phases.npz 256 wnoise_music_256.bin   # level 8
# edit ref_center / ref_extent for your target region, then
MUSIC examples/music_zoom.conf           # → ic_quaia_zoom.gdt
```
This is the cleanest setup: our reconstruction supplies the constrained large scales, MUSIC supplies the
zoom — both use the same multi-scale white-noise construction, so they compose exactly.

## 4. N-GenIC / Gadget / 2LPTic

These codes generate their *own* phases from a seed and do not read an external white-noise field, so
there is no direct drop-in. Two supported routes:

1. **Via MUSIC** — produce the Gadget IC with MUSIC (§3, `format = gadget2`), which N-body codes read
   directly. This is the recommended path.
2. **Via DiscoDJ.jl** — generate the displacement/particle IC yourself (§5) and write your code's format.

(If you maintain a patched N-GenIC that ingests external white noise, the field layout is documented in
`examples/make_music_wnoise.jl`; match the unit-variance, mean-zero, column-major convention.)

## 5. DiscoDJ.jl directly (pure Julia, no external IC code)

```bash
julia examples/snapshot_discodj.jl quaia_icbox_phases.npz 512 49.0    # res, z_init
```
This applies the IC operator + 2LPT at z_init and returns periodic particle positions — the same forward
DiscoDJ uses, so the constrained field is preserved exactly. The example writes positions; add 2LPT
velocities and your header for a ready-to-run snapshot, or use DiscoDJ's writers.

## 6. Conventions (match these when ingesting the field)

- **Units:** comoving Mpc/h. **Cosmology:** Planck18-like (`manifest.json`).
- **Grid:** periodic, `res^3`, column-major (x fastest, then y, then z). **White noise:** unit variance,
  zero mean, fixed amplitude (|ω̂(k)| = const).
- **Gauge:** −1/k² (infall into overdensities). **Observer:** box centre.
- **Sign:** verify on first use (MUSIC negates by default; `grafic_sign = yes` disables that). The phase
  convention is correct iff overdensities land where the Quaia quasars are.

Questions / issues: open an issue on this repo or on
[DiscoInverse.jl](https://github.com/yipihey/DiscoInverse.jl).
