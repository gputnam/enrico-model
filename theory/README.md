# nuflow theory notes

Each physics component in `nuflow` has a note here describing the model, the equations it implements, the approximations made, and how it is validated. The notes use GitHub-flavoured Markdown with inline LaTeX.

| Note | Component(s) |
|---|---|
| [framework.md](framework.md) | Staged sampling algorithm, cross-section tables, unweighting, normalization, reproducibility |
| [quasielastic.md](quasielastic.md) | Stage 1 `CCQEChannel`, stage 2 `RelativisticFermiGas` / `FreeNucleon`, stage 3 `QuasiElasticKinematics`, `FormFactors`, `QECurrent` |
| [spectral_function.md](spectral_function.md) | Stage 2 `ShellModelSpectralFunction` (12C HO shells + (e,e'p) removal energies + SRC np pairs with partner emission, local Pauli blocking, vertex position handed to the cascade) |
| [local_fermi_gas.md](local_fermi_gas.md) | Stage 2 `LocalFermiGas` (Q-value energy balance) and stage 3 `ValenciaRPA` (Nieves–Amaro–Valverde RPA, `rpa: true`), `apps/nuflow-lfg-validate` |
| [pion_production.md](pion_production.md) | Stage 1 `CCRESChannel`, stage 3 `ResonantPionKinematics`, `PionProductionCurrent` (HNV Δ(1232) + chiral background), `core/Dirac.h` |
| [nstar.md](nstar.md) | Stage 1 `CCNStarChannel`, stage 3 `NStarKinematics` (N(1440), N(1520), N(1535): Lalakulich–Paschos–Piranishvili form factors, energy-dependent Breit–Wigners, N π / N η / N π π decays) |
| [dis.md](dis.md) | Stage 1 `CCDISChannel`, stage 3 `DISKinematics`, `BodekYang` (effective-LO structure functions), `GRV98LO`, `DISHadronization` |
| [mec.md](mec.md) | Stage 1 `CCMECChannel`, stage 3 `MECKinematics`, `MECCurrent` (relativistic 2p2h MEC of Ruiz Simo *et al.*: seagull, pion-in-flight, pion-pole, Δ), `apps/nuflow-mec-response` |
| [fsi_none.md](fsi_none.md) | Stage 4 `TransparentNucleus` |
| [fsi_cascade.md](fsi_cascade.md) | Stage 4 `IntranuclearCascade` (nucleon + pion INC: local density, Cugnon NN, Oset–Salcedo Δ-hole pion absorption, SKAT formation zone), `NuclearDensity`, `HadronNucleonXSec`, `apps/nuflow-fsi-validate` |

## Conventions

- **Units:** natural units, $\hbar=c=1$. Energies and momenta in GeV; cross sections in $10^{-38}\,\text{cm}^2$.
- **Metric and Dirac algebra:**
  - metric $g=\mathrm{diag}(+,-,-,-)$;
  - $\gamma_5=i\gamma^0\gamma^1\gamma^2\gamma^3$;
  - Levi-Civita symbol with $\epsilon^{0123}=-1$, used only in the derivations.
- **Kinematics:**
  - $k$ and $k'$ are the neutrino and charged-lepton four-momenta; $p$ and $p'$ the initial and final nucleon four-momenta;
  - $q=k-k'$ and $Q^2=-q^2$.
- **PDG codes:** particles and nuclei use PDG Monte Carlo codes (nuclei as `10LZZZAAAI`).
- **Particle status codes:** 4 beam, 20 target, 21 struck nucleon, 22 hadron before FSI, 23 nuclear remnant (before FSI), 25 nuclear remnant after the cascade, 1 final state.

## Data provenance

| File | Source |
|---|---|
| `data/flux/minerva_le_2016/` | L. Aliaga *et al.* (MINERvA), *Phys. Rev. D* **94**, 092005 (2016) [arXiv:1607.00704], ancillary files |
| `data/flux/minerva_me_2019/` | E. Valencia *et al.* (MINERvA), *Phys. Rev. D* **100**, 092001 (2019) [arXiv:1906.00111], ancillary files |
| `data/pdf/grv98lo.dat` | GRV98 LO grid: M. Glück, E. Reya, A. Vogt, *Eur. Phys. J. C* **5**, 461 (1998) [hep-ph/9806404] (see [dis.md](dis.md)) |
| `data/electron/qes_12C.dat` | 12C inclusive $(e,e')$ cross sections from the quasielastic electron-nucleus scattering archive (O. Benhar, D. Day, I. Sick, *Rev. Mod. Phys.* **80**, 189 (2008); arXiv:nucl-ex/0603032), used for validation only |
| `data/nuclei/nuclei.yaml` | AME2020 masses; CC Q values (Nieves *et al.* 2004, Table I); Fermi-gas parameters from E. J. Moniz *et al.*, *Phys. Rev. Lett.* **26**, 445 (1971), and common generator choices (see the file header) |
