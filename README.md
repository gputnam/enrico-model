# nuflow

`nuflow` is a staged event generator for GeV-scale neutrino–nucleus interactions. Every event passes through four stages, and each stage is an abstract C++ class that can be swapped for another implementation:

| Stage | Directory | Base class | Implementations (v0.1) |
|---|---|---|---|
| 1. Channel and target selection | `channel/` | `InteractionChannel` (+ `ChannelSelector`) | `CCQEChannel`, `CCRESChannel`, `CCNStarChannel`, `CCDISChannel`, `CCMECChannel`, `CCCOHChannel` |
| 2. Initial-state nuclear sampling | `initial_state/` | `NuclearModel` | `RelativisticFermiGas`, `LocalFermiGas`, `ShellModelSpectralFunction` (12C SF + SRC pairs), `FreeNucleon`, `GroundStateNucleus` (coherent) |
| 3. Kinematics generation | `kinematics/` | `KinematicsGenerator` | `QuasiElasticKinematics` (+ `FormFactors`), `ResonantPionKinematics` (+ `PionProductionCurrent`), `DISKinematics` (+ `BodekYang`, `GRV98LO`, `DISHadronization`), `MECKinematics` (+ `MECCurrent`), `CoherentPionKinematics` (+ `PionNucleusElastic`) |
| 4. Final-state interactions | `fsi/` | `FSIModel` | `TransparentNucleus` (no FSI), `IntranuclearCascade` (+ `NuclearDensity`, `HadronNucleonXSec`) |

The `orchestration/` layer does the rest:
- builds the total-cross-section tables, caching them on disk;
- samples $E_\nu$ from $\Phi(E)\,\sigma_\text{tot}(E)$;
- runs the four stages on worker threads;
- writes the output.

Shared types (four-vectors, RNG, event record, nuclei, configuration) live in `core/`, and the flux sources and readers in `flux/`.

The physics behind each component is written up in [`theory/`](theory/README.md). The quasi-elastic model is described in [`theory/quasielastic.md`](theory/quasielastic.md), single-pion production in the Δ(1232) region (Hernández–Nieves–Valverde Δ + chiral background, process id 400) in [`theory/pion_production.md`](theory/pion_production.md), the second-resonance-region states $N(1440)$, $N(1520)$, $N(1535)$ (Lalakulich–Paschos–Piranishvili form factors, process id 401) in [`theory/nstar.md`](theory/nstar.md), and inclusive shallow/deep-inelastic scattering above $W=1.6$ GeV (Bodek–Yang effective-LO structure functions on GRV98 LO, hadronized into a nucleon and pions, process id 600) in [`theory/dis.md`](theory/dis.md), and 2p2h excitations by relativistic meson-exchange currents (Ruiz Simo *et al.* 2017, two nucleons out, process id 300) in [`theory/mec.md`](theory/mec.md), and charged-current coherent pion production on nuclei (Berger–Sehgal PCAC model with the lepton-mass correction and π–¹²C elastic data, process id 100) in [`theory/coherent_pion.md`](theory/coherent_pion.md). Final-state interactions are a semi-classical intranuclear cascade for nucleons and pions (local-density ¹²C, Cugnon NN cross sections, Oset–Salcedo Δ-hole pion absorption, SKAT formation zone) in [`theory/fsi_cascade.md`](theory/fsi_cascade.md).

## Build

Everything comes from the spack stack on cvmfs:
- **Required:** ROOT 6.28, yaml-cpp, gcc 12.
- **Optional:** HepMC3, for NuHepMC output.
- **Tests only:** Catch2.

```bash
source setup.sh                       # spack load of the pinned toolchain + NUFLOW_DIR
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j16
./build/nuflow-tests                  # unit + physics validation tests
```

## Run

```bash
# 1M events, MINERvA LE FHC numu flux on CH (default config: 16 threads)
./build/nuflow-gen --config config/minerva_le_fhc_numu_CH.yaml -o runs/minerva_le_numu_CH.root

# mono-energetic beam, other species/targets, NuHepMC output
./build/nuflow-gen --config config/mono_1GeV_C12.yaml --mono 2.0 --nu -14 --target 1000080160:1 \
                   -n 100000 --hepmc out.hepmc3

# print the tabulated total cross sections
./build/nuflow-xsec --config config/mono_1GeV_C12.yaml --emin 0.2 --emax 20 --nodes 40

# validation plots
root -l -b -q 'scripts/validate.C("runs/minerva_le_numu_CH.root", "runs/plots")'
```

### Configuration

Runs are configured in YAML (see `config/`).

| Section | Contents |
|---|---|
| `run` | Events, seed, threads, output, `batch_size`, `safety_factor`, ROOT `compression` |
| `beam.pdg` | Neutrino species |
| `beam.flux` | `{type: mono, energy: E}`, `{type: flat, emin, emax}`, or a histogram (see below) |
| `target` | A list of `{pdg, count}` or the string `"1000060120:1,1000010010:1"` |
| `channels` | A list of channels with their stage-2 and stage-3 configuration (see below) |
| `fsi` | The stage-4 model: `{type: none}` (transparent) or `{type: cascade, ...}` (see below) |
| `xsec_table` | Energy nodes, MC samples per node, cache directory |

Channel types:
- `CCQE`: `nuclear_model` (`RFG`, `LFG` ([theory/local_fermi_gas.md](theory/local_fermi_gas.md)), `SF` (12C spectral function with SRC partner emission, [theory/spectral_function.md](theory/spectral_function.md))), `form_factors: {vector: BBBA05|dipole, axial: dipole|zexp, MA, gA, ...}`, `rpa` (Valencia RPA, `true` or a parameter map; default off), `proposal_q2_scale`. `LFG`/`SF` are CCQE-only; the other channels use `RFG`.
- `CCRES` (single pion, Δ region, [theory/pion_production.md](theory/pion_production.md)): `nuclear_model`,
  `delta: {C5A0: 1.00, MA: 0.93, mass: 1.232, f_star: 2.14}`, `background: true|false`,
  `nucleon_form_factors` (background terms; default dipole/Galster, gA 1.26, MA 1.05), `w_max` (default 1.6 GeV),
  `proposal_q2_scale`, `proposal_bw_fraction`.
- `CCNSTAR` (N(1440), N(1520), N(1535) for nu n / nubar p, [theory/nstar.md](theory/nstar.md), process id 401):
  `nuclear_model`, `w_max` (default 1.6 GeV, = CCRES `w_max`), `MV` (0.84), `MA` (1.05),
  `resonances` (default `[P11_1440, D13_1520, S11_1535]`), `proposal_q2_scale`.
- `CCDIS` (inclusive inelastic above `w_min`, [theory/dis.md](theory/dis.md)): `nuclear_model`, `w_min` (default 1.6 GeV;
  must equal the CCRES `w_max`), `w_max` (default 0 = kinematic limit), `bodek_yang: {A: 0.621, B: 0.380, N: 1.026,
  axial: typeII|typeI, charm: true, du_correction: true}`, `hadronization: {dispersion_slope: 0.35, mean_pt2: 0.25,
  w_lps: 2.3}`, `single_pion: {scale_i32, scale_i12, isospin}` (isospin-decomposed rescaling of the single-pion final
  states, inclusive unchanged; dis.md section 3.1; default none), `proposal_q2_scale`.
  Needs the PDF grid `data/pdf/grv98lo.dat`.
- `CCMEC` (2p2h meson-exchange currents, [theory/mec.md](theory/mec.md)): `nuclear_model`,
  `mec: {C5A0: 1.00, MA: 0.93, C3V0: 2.13, MV: 0.84, f_star: 2.13, delta_propagator: real|full, lambda_piNN: 1.3,
  lambda_piNDelta: 1.15, mass: 1.232, gA: 1.26, f_pi: 0.093}`, `nucleon_form_factors` (F1V only), `proposal_q_scale`
  (default 1 GeV), `q3_max` (default 0 = no cut).
- `CCCOH` (CC coherent pi+- on nuclei, A >= 4, [theory/coherent_pion.md](theory/coherent_pion.md), process id 100; no
  `nuclear_model`): `f_pi` (0.1307), `mA` (0.95), `lepton_mass_correction` (true),
  `elastic: {R0: 1.057, table_max_T: 0.87, rs_min_T: 1.7}` (pi-12C data below `table_max_T`, Rein-Sehgal hadronic
  model above `rs_min_T`), `proposal_slope` (0.8). The pion is written as a final-state particle (not cascaded) and the
  recoiling ground-state nucleus as the remnant.

FSI: `{type: cascade}` ([theory/fsi_cascade.md](theory/fsi_cascade.md)) options `step` (0.2 fm), `nucleon_rescattering`, `nn_pion_production`,
`pion_quasielastic`, `pion_absorption` (all true), `like_pair_absorption` (0.05), `formation_zone` (true), `formation_mu2` (0.08 GeV²),
`formation_processes` ([600]). Needs `rch`, `Sp`, `Sn` in `data/nuclei/nuclei.yaml`.

A histogram flux is given as `{type: histogram, format: root|csv|minerva_le_txt|table, file, histogram, density: per_bin|per_GeV, scale, units, emin, emax}`.

Flux formats that are supported:
- **ROOT `TH1`:** any histogram, selected by name. Example: the MINERvA LE 2016 release, `data/flux/minerva_le_2016/minerva_flux.root`, histograms `numu_fhc` etc. They are in ν/m²/10⁶ POT per 0.5 GeV bin, so use `density: per_bin, scale: 1e-6`.
- **MINERvA CSV:** `BinLow,BinHigh,Flux`, as in the ν–e-constrained ME 2019 release in `data/flux/minerva_me_2019/` (ν/m²/POT/GeV).
- **MINERvA LE text:** the same release as the ROOT file, as `*.txt`.
- **Generic table:** `E_low E_high flux` or `E_center flux`.

### Output

The ROOT file contains:
- **`events` TTree**, one entry per event. Event-level branches are `Enu`, `Q2`, `q0`, `q3`, `W`, `channel`, `target_pdg`, the struck nucleon (`hit_*`) and `weight`. The particle arrays are `pdg`, `status`, `mother`, `E`, `px`, `py` and `pz` (up to 256 particles).
- **`run` TTree** with the flux-averaged cross section, trial and overweight counters, and the seed.
- **σ(E) graphs** per channel and nucleus, the flux histogram, and the full configuration.

The particle status codes follow NuHepMC:

| Status | Meaning |
|---|---|
| 4 | Beam |
| 20 | Target |
| 21 | Struck nucleon |
| 22 | Hadron before FSI |
| 23 | Nuclear remnant (before FSI) |
| 25 | Nuclear remnant after the cascade (excited; carries the balancing four-momentum) |
| 1 | Final state |

`--hepmc` additionally writes a NuHepMC v1.0 HepMC3 ASCII file.

Events are **unweighted** (`weight = 1`). The absolute normalization is set by `flux_avg_xsec`, in 10⁻³⁸ cm² per target unit (e.g. per CH):

$$\frac{d\sigma}{dX} = \texttt{flux\_avg\_xsec}\times\frac{N(X)}{N_\text{events}\,\Delta X}.$$

### Reproducibility and performance

- **Reproducibility:** events are generated in batches. Batch $b$ uses RNG stream $(\text{seed}, b)$ and batches are written in order, so the output depends only on the seed, never on the thread count.
- **Speed:** a trial of stages 2 and 3 takes about 0.2 µs, and the QE acceptance is about 50%.
  - Without I/O, generation runs at 1.3×10⁶ events/s on one thread and 6.6×10⁶ events/s on 16 threads.
  - With output, throughput is set by the ROOT writer: about 2.8×10⁵ events/s. The 1M-event MINERvA test takes 3.6–4.6 s of generation, depending on machine load (6–7 s wall, including building the cross-section table).

## Repository layout

```
core/            event record, four-vectors, Dirac matrices, RNG, PDG helpers, nuclei, YAML config
flux/            FluxSource, MonoEnergeticFlux, HistogramFlux, flux file readers
channel/         stage 1: InteractionChannel, CCQEChannel, CCRESChannel, CCNStarChannel, CCDISChannel, CCMECChannel, CCCOHChannel,
                 ChannelSelector
initial_state/   stage 2: NuclearModel, RelativisticFermiGas, LocalFermiGas, ShellModelSpectralFunction, FreeNucleon,
                 GroundStateNucleus
kinematics/      stage 3: KinematicsGenerator, QuasiElasticKinematics, FormFactors, QECurrent, ResonantPionKinematics,
                 PionProductionCurrent, DISKinematics, BodekYang, GRV98LO, DISHadronization, MECKinematics, MECCurrent,
                 CoherentPionKinematics (+ PionNucleusElastic)
fsi/             stage 4: FSIModel, TransparentNucleus, IntranuclearCascade, NuclearDensity, HadronNucleonXSec
orchestration/   Generator, XSecTableBuilder, EnergySampler, ModelFactory, writers
apps/            nuflow-gen, nuflow-xsec, nuflow-qe-table, nuflow-mec-response (2p2h response functions),
                 nuflow-fsi-validate (hadron-nucleus and event-level FSI validation)
config/          example run configurations
data/            flux files (with provenance READMEs), nuclear parameters, GRV98 LO PDF grid (data/pdf), xsec cache
theory/          physics documentation per model
tests/           Catch2 unit and physics-validation tests
scripts/         validation plotting, run_loop_prediction.sh (closed loop)
```

## Adding a model

1. Derive from the base class of the stage you are extending.
   - For a channel, implement `Applies`, `Request` (which nucleons and how many), `Multiplicity` and `ProcessId`. A channel acting on either nucleon species overrides `SampleRequest` to draw the species per trial (see `CCRESChannel`).
   - For a nuclear model, implement `Sample` and `PauliBlocked`.
   - For a kinematics generator, implement `Generate`. It returns the weight $d\sigma/(\text{proposal density})$ in 10⁻³⁸ cm² per struck nucleon (or pair). Weight-independent completion of the final state (e.g. DIS hadronization) goes in `Finalize`, which runs only for accepted events.
   - For an FSI model, implement `Propagate`.
2. Register the new type in `orchestration/ModelFactory.cc`.
3. Add a theory note in `theory/` (and index it in `theory/README.md`).
4. Add the new `.cc` files to the `add_library` list in `CMakeLists.txt`, and a test in `tests/`.
5. Use a NuHepMC-conformant `ProcessId` (100 coherent, 200 CCQE, 300 CC MEC/2p2h, 400 CC RES,
   500 CC SIS/non-resonant, 600 CC DIS; NC = +50): the measurement bridge maps it to the GENIE
   interaction type (`measurement/xsec/xsec/genpred.py`).

The cross-section tables and the unweighting are generic: they come from the channel's own trial weight, so a new channel needs no orchestration changes.

**Cross-section cache.** `data/cache/` is keyed on the channel's YAML node. A *code* change that
alters the trial weight for an unchanged config is not detected: bump `kCacheVersion` in
`orchestration/XSecTableBuilder.cc` (or `rm -rf data/cache/*`).

## Closed-loop usage (orchestrate/)

The closed loop compares this model with the MINERvA measurements in `measurement/`:

```bash
scripts/run_loop_prediction.sh <workspace>/model.root     # build + tests + 2M events, ~15 s
```

uses `config/minerva_me_fhc_numu_CH_loop.yaml`: νμ on CH with the *measurement* flux
(`data/flux/minerva_me_fhc_measurement/`, the data-POT-weighted ν–e constrained ME FHC flux the
measurements divide by).  Change the model in the `channels`/`fsi` sections (and in code); keep the
`run`, `beam` and `target` sections fixed so predictions remain comparable between iterations.
Predictions are normalised per nucleon (`flux_avg_xsec / nucleons_per_target`) and evaluated in each
measurement's own signal definition, so every final-state particle written with status 1 matters:
the measurements classify events by their final-state topology (muon, pions, protons, photons).
