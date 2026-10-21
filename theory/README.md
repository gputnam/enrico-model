# nuflow theory notes

Each physics component in `nuflow` has a note here describing the model, the equations it implements, the approximations made, and how it is validated. The notes use GitHub-flavoured Markdown with inline LaTeX.

| Note | Component(s) |
|---|---|
| [framework.md](framework.md) | Staged sampling algorithm, cross-section tables, unweighting, normalization, reproducibility |
| [quasielastic.md](quasielastic.md) | Stage 1 `CCQEChannel`, stage 2 `RelativisticFermiGas` / `FreeNucleon`, stage 3 `QuasiElasticKinematics`, `FormFactors`, `QECurrent` |
| [pion_production.md](pion_production.md) | Stage 1 `CCRESChannel`, stage 3 `ResonantPionKinematics`, `PionProductionCurrent` (HNV Δ(1232) + chiral background), `core/Dirac.h` |
| [dis.md](dis.md) | Stage 1 `CCDISChannel`, stage 3 `DISKinematics`, `BodekYang` (effective-LO structure functions), `GRV98LO`, `DISHadronization` |
| [mec.md](mec.md) | Stage 1 `CCMECChannel`, stage 3 `MECKinematics`, `MECCurrent` (relativistic 2p2h MEC of Ruiz Simo *et al.*: seagull, pion-in-flight, pion-pole, Δ), `apps/nuflow-mec-response` |
| [fsi_none.md](fsi_none.md) | Stage 4 `TransparentNucleus` |

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
- **Particle status codes:** 4 beam, 20 target, 21 struck nucleon, 22 hadron before FSI, 23 nuclear remnant, 1 final state.

## Data provenance

| File | Source |
|---|---|
| `data/flux/minerva_le_2016/` | L. Aliaga *et al.* (MINERvA), *Phys. Rev. D* **94**, 092005 (2016) [arXiv:1607.00704], ancillary files |
| `data/flux/minerva_me_2019/` | E. Valencia *et al.* (MINERvA), *Phys. Rev. D* **100**, 092001 (2019) [arXiv:1906.00111], ancillary files |
| `data/pdf/grv98lo.dat` | GRV98 LO grid: M. Glück, E. Reya, A. Vogt, *Eur. Phys. J. C* **5**, 461 (1998) [hep-ph/9806404] (see [dis.md](dis.md)) |
| `data/nuclei/nuclei.yaml` | AME2020 masses; Fermi-gas parameters from E. J. Moniz *et al.*, *Phys. Rev. Lett.* **26**, 445 (1971), and common generator choices (see the file header) |
