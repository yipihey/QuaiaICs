# Method — how the Quaia-constrained ICs are derived

This note gives the full derivation behind the boxes in this repository. The implementation lives in
[DiscoInverse.jl](https://github.com/yipihey/DiscoInverse.jl) (`reconstruct_quaia`, `constrained_ic_box`,
`refine_phases`) on top of the differentiable ΛCDM engine
[DiscoDJ.jl](https://github.com/yipihey/DISCO-DJ).

## 1. The data

[Quaia](https://arxiv.org/abs/2306.17749) (Storey-Fisher et al. 2024) is an all-sky Gaia DR3 × unWISE
quasar catalog. We use the G < 20.0 sample with the Galactic-plane cut (|b| ≥ 10°), 744,834 quasars over
z ≈ 0–4.7. Each quasar has an accurate sky position (ra, dec) and a *spectro-photometric* redshift with
error σ_z ≈ 0.03–0.12. The radial blur is

  σ_χ = (∂χ/∂z) σ_z ≈ 120 Mpc/h  (median),

orders of magnitude larger than the clustering scale — so the line-of-sight map is badly smeared while
the transverse map is sharp.

## 2. The differentiable forward model

A primordial unit-variance white-noise field `ω` on a periodic grid is mapped to a model quasar density
through the standard ΛCDM chain, every step differentiable (Zygote), implemented in DiscoDJ.jl:

  ω  →  φ(k) = −√(P(k))/k² · ω̂(k)            [IC operator: fixed cosmology transfer]
     →  ψ(q) = nLPT displacements from φ        [1–3 LPT, growth factors]
     →  x_obs(q) on the past lightcone          [observer at box centre, a(χ)]
     →  ρ_g(x) tetrahedral CDM-sheet density    [bias-weighted, grid-free at the tracer points]

The quasar likelihood is an inhomogeneous Poisson point process in `ρ_g`, restricted to the survey
window, with a linear bias b₁ = 2.5.

## 3. Joint inference of the field and the true distances

We **parametrize each quasar by its comoving radial distance** χ_i rather than redshift, so the embedding
`x_i = χ_i n̂_i + shift` is linear in the free parameter (differentiating z→χ through `comoving_distance`
is numerically unstable). The quasar positions become free, differentiable parameters via a query-point
density gradient ∂ρ_g/∂x = ∇ρ. We then minimize

  L(ω, {χ_i}) = −Σ_i log ρ_g(x_i)  +  Σ_i (χ_obs,i − χ_i)² / (2 σ_χ,i²)  +  (field prior),

i.e. the 3-D clustering pulls each quasar's noisy radial position toward the cosmic web, regularized by
its photo-z. The reconstruction tightens the per-object radial uncertainty from ~120 Mpc/h to a few Mpc/h
where the field is informative, and leaves it photo-z-wide in voids and gaps (figure in the README).
RSD (~2 Mpc/h) ≪ σ_χ, so the forward is run in real space.

## 4. Fixed-amplitude phases → a constrained realization

The field is parametrized by its **Fourier phases only**:

  ω(φ) = irfft( e^{iφ} ) / std,    |ω̂(k)| = const,

the Angulo–Pontzen *fixed-amplitude* construction. Every mode amplitude is pinned to √P(k) through the
IC operator, so the field has **exactly** the cosmological P(k) — no per-mode amplitude scatter. The
inference (an alternating MAP over φ and {χ_i}) moves the phases the Quaia data constrain; the rest keep
their random draw. The output is therefore a **constrained realization**: a legitimate random ΛCDM field
that *also* reproduces the Quaia clustering within the photo-z errors. Different random seeds for the
unconstrained phases form an ensemble (same constrained structure, independent elsewhere).

This is an MAP estimate of the constrained modes, not a posterior sample — the honest uncertainty is the
ensemble spread over seeds.

## 5. Coarse-constrain → fine-realize (arbitrary resolution)

The differentiable forward is memory-bound and caps near 384³ on a single GPU. But §1 shows the Quaia
constraint carries **no small-scale information** — the photo-z blur band-limits it to k ≲ 0.05 h/Mpc.
So we split the resolution:

**Constrain (coarse).** Run the full reconstruction on a grid (res = 256, dx ≈ 52 Mpc/h) whose Nyquist
comfortably exceeds the informed scale. This loses nothing and produces the constrained phases φ₂₅₆ — the
compact **carrier** (`quaia_icbox_phases.npz`, ~67 MB).

**Realize (fine).** Build the box at any resolution N by **spectral white-noise refinement**:

  1. Allocate fresh random phases on the fine rfft grid (N/2+1, N, N).
  2. *Embed* φ₂₅₆ into it at matching physical wavenumber — for a uniform-L refinement the coarse mode at
     integer frequency n maps to the same n on the fine grid (`refine_phases` / `_embed_indices`, which
     handle the rfft half-axis and the Nyquist/Hermitian wrapping). Interior modes transfer bit-exactly.
  3. ω_N = phase_field(φ_N): unit modulus everywhere ⇒ exactly fiducial P(k) at all scales, the
     constrained modes carried at large scales, fresh uniform-random phases filling the small scales.

This is a pure FFT — no forward model, no autodiff tape — so it scales to 1024³, 2048³, … limited only by
ordinary array memory. It is the same multi-scale white-noise idea as MUSIC (Hahn & Abel 2011), here
seeded by a field-level constraint.

**Validation (measured on the shipped box):**

- decoupling: a 1024³ box reproduces the res-256 constraint with cross-correlation **1.0000**;
- fixed P(k): interior per-mode |ω̂(k)|² constant to **std/mean ≈ 9×10⁻¹⁶** (only the kx=0 / Nyquist
  reality planes deviate, as for any real fixed-amplitude field);
- ensemble: two realizations share every mode below the coarse Nyquist (r = 1) and decorrelate above it
  (r = 0) — see `figures/powerspectrum.png`;
- a 2048³ box (140× more cells than the forward can hold) realizes in ~7 min on CPU.

## 6. Parameters for this carrier

| parameter | value |
|---|---|
| cosmology | Ω_m=0.315, Ω_b=0.049, h=0.674, σ₈=0.81, n_s=0.965 |
| box size L | 13 260 Mpc/h (= 1.3 × the all-sky Quaia extent) |
| res_constrain | 256 (dx ≈ 52 Mpc/h) |
| observer | box centre, (6640, 6649, 6620) Mpc/h |
| constrained sphere | radius χ(z_max=4.66) ≈ 5 230 Mpc/h |
| tracers | 744,834 Quaia quasars, b₁ = 2.5, real space |
| fixed amplitude | yes (Angulo–Pontzen) |
| seed (carrier) | 101 |

See `manifest.json` for the machine-readable version, and
[DiscoInverse.jl](https://github.com/yipihey/DiscoInverse.jl) for the source.
