# Data — the constraint carrier

The data product is a single compact file:

### `quaia_icbox_phases.npz`  (~67 MB)

The **constraint carrier** — the Quaia-constrained fixed-amplitude phases on the res = 256 grid
(`phi_coarse`, shape `(129, 256, 256)`, complex/real per NPZ) plus the full manifest. **It fully
determines the box**: every resolution (1024³, 2048³, …) is re-realized from it locally with
`examples/realize_box.jl`, so we do not ship multi-GB grids.

**It is committed in this repository** at `data/quaia_icbox_phases.npz`, so a plain `git clone` gives you
everything. (If you prefer to keep it out of your working tree, it can also be served as a release
asset — `gh release download ... --pattern 'quaia_icbox_phases.npz'`.)

Then, e.g.:

```bash
julia examples/realize_box.jl data/quaia_icbox_phases.npz 1024 omega_1024.f32
```

### Contents

| key | meaning |
|---|---|
| `phi_coarse` | constrained fixed-amplitude phases, rfft grid `(res/2+1, res, res)` at res = 256 |
| `boxsize`, `res_constrain`, `res_box` | geometry (13260 Mpc/h, 256, 1024) |
| `Omega_m, Omega_b, h, sigma8, n_s` | cosmology |
| `observer`, `shift` | observer / comoving→box offset (box centre) |
| `constrained_radius` | χ(z_max) ≈ 5229 Mpc/h — the constrained sphere |
| `z_min, z_max, n_quasars, b1, seed, fixed_amplitude` | constraint provenance |

See [`../manifest.json`](../manifest.json) for the machine-readable summary and
[`../docs/method.md`](../docs/method.md) for how it was produced.

> Realized white-noise grids and particle snapshots are **not** committed (a 1024³ float32 box is 4 GB,
> a 2048³ is 34 GB). Generate them locally from the carrier with the example scripts.
