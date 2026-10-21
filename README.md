# PUFB–OTFS (R12) — MATLAB Simulation Scripts

Code-only package with the **31 MATLAB scripts** that produce the numerical results
reported in the *PUFB × OTFS* manuscript (revision R12). No data, figures, LaTeX,
or Python files are included.

---

## 1. Requirements

- MATLAB **R2024a** or later (developed and validated on R2024a Update 4)
- **Signal Processing Toolbox**
- **Statistics and Machine Learning Toolbox**
- *Parallel Computing Toolbox* (optional — accelerates `parfor` loops; without it
  they run serially)

## 2. Package layout (flat, numbered by pipeline order)

```
PUFB-OTFS-R12-MATLAB-scripts/
├── README.md
├── 01_channel_generator/        generate_channel_r20a.m     Physical 2D delay–Doppler channels (EVA/ABID)
├── 02_shared_lib/               fir_mallat_r22a.m, r22a_io.m  Mallat FIR bank + I/O & SHA-256 helpers
├── 03_axis_specific_pufb/       run_r22d.m + r22d_*.m       R22D  — receiver-aware axis-specific PUFB
├── 04_ber_confirmation/         run_r22dcb.m + r22dcb_*.m   R22DCB — BER confirmation (phases A/B)
├── 05_highsnr_mismatch_audit/   run_r22dhma.m + r22dhma_*.m R22DHMA — high-SNR mismatch diagnostic
├── 06_receiver_floor/           run_r22e.m + r22e_*.m       R22E  — normalized-loading receiver floor
├── 07_natural_baselines/        run_r22h.m                  R22H  — natural baseline transforms
└── 08_main_campaign/            run_r29.m + r29_*.m         R2.9  — main PUFB design + OTFS campaign
```

Each `run_*.m` is the entry point of its module; the other files in the same folder
are its helpers.

## 3. Quick start

```matlab
% 1) Put ALL package code on the path (one line):
cd PUFB-OTFS-R12-MATLAB-scripts
addpath(genpath(pwd));

% 2) Tell the scripts where the project data live:
root = '/path/to/repository';   % folder containing pufb_otfs_v2/results/...

% 3) Run any stage of any module:
run_r29(root,'CONTROLS');       % fast sanity controls (paraunitarity, orthogonality)
run_r22d(root,'PREPARE');       % then TRAIN -> FREEZE -> TEST, in this order
```

All runners have the same signature: `result = run_XXX(root, stage)`.

## 4. How paths work (code vs. data)

The scripts locate **their own code** automatically (via `mfilename`), so this
package can live anywhere and be renamed freely. The argument `root` is used only
for **data, configs and results**, which stay in the original repository layout:

- read/written under `root/pufb_otfs_v2/results/...` (seed plans, manifests,
  channel banks, output tables);
- preregistration/config JSONs read from `root/pufb_otfs_v2/audit/...`;
- four helper directories that are *not* part of this package are also taken from
  the repository: `r2_2c_confirm_strict`, `r2_2c_receiver_aware_fixed_budget`,
  `r2_5_abid_full_mmse`, `r2_2g_factorial_ablation`, `r1_0r_wotfs_reproduction`.

So: **code from this package, data from the repository.** Point `root` at the
repository (or at any folder that reproduces its `pufb_otfs_v2/results` tree).

## 5. Run order and stage gates

| Module | Entry point | Stages (in order) | Depends on |
|---|---|---|---|
| R22D | `run_r22d.m` | `PREPARE → TRAIN → FREEZE → TEST` | — |
| R22DCB | `run_r22dcb.m` | `PREPARE → A → B` | frozen R22D outputs |
| R22DHMA | `run_r22dhma.m` | `PREPARE → RUN` | frozen `r2_2d_ber_ext` outputs |
| R22E | `run_r22e.m` | `PREPARE → DEVELOP → FREEZE → TEST → BER` | frozen R22DCB outputs |
| R22H | `run_r22h.m` | `PREPARE → TEST` | frozen R22D + R22G outputs |
| R2.9 | `run_r29.m` | `CONTROLS → DIMAUDIT → DATA → DESIGN → REFINE → VALIDATE → TEST → STATS → BER` | seed plan in results tree |

Calling a stage out of order is blocked by built-in temporal-gate assertions:
TEST channels cannot be generated or accessed before the corresponding
candidate/support freeze manifests exist.

## 6. Reproducibility safeguards

- Deterministic seed plans (physical channels, data, noise and bootstrap use
  independent documented bases);
- SHA-256 verification of code and input files before/after every stage;
- Channel-/SNR-indexed randomness — interrupted runs resume bit-exactly.

> **Note on hashes:** the directories of this package were renamed and the path
> lines of the six `run_*.m` files were adapted accordingly. The scientific code
> is byte-identical to the frozen version, but the `run_*.m` files themselves no
> longer match the frozen code-hash manifests. The scripts recompute and verify
> hashes at runtime; expect the recorded `run_*.m` hashes to differ from the
> original manifests.

## 7. Notes

- The few Code Analyzer messages in the `r29_*` files are benign style/performance
  suggestions inherited from the production scripts; the scientific code was left
  untouched.
- Tested on macOS (Apple Silicon), MATLAB R2024a Update 4. On macOS, launch with
  the native architecture: `matlab -maca64`.
