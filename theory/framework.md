# Framework: staged event generation and its normalization

This note explains why the four-stage algorithm in `orchestration/Generator.cc` produces events with the correct distribution, and how the output is normalized.

## 1. What is generated

Notation:
- The target unit (e.g. one CH "molecule") contains $`n_A`$ nuclei of each species $`A`$.
- The neutrino flux has spectral density $`\Phi(E)`$; only its shape matters for event generation.
- Each channel $`c`$ (e.g. CCQE) acting on nucleus $`A`$ has the fully differential cross section
```math
d\sigma_{c,A}(E) = \int d\mathcal{I}_A \; \frac{d\sigma_c}{d\mathcal{K}}(E,\mathcal{I})\, d\mathcal{K},
```
where:
- $`\mathcal{I}`$ is the nuclear initial state (struck-nucleon momenta, removal energies);
- $`d\mathcal{I}_A`$ is its probability measure under the nuclear model, times the number of struck-nucleon configurations;
- $`\mathcal{K}`$ are the final-state kinematic variables.

The density of events the generator must produce is
```math
p(E, A, c, \mathcal{I}, \mathcal{K}) \;\propto\; \Phi(E)\; n_A\; \frac{d\sigma_c}{d\mathcal{K}}(E, \mathcal{I})\; d\mathcal{I}_A\, d\mathcal{K} ,
\qquad(1)
```
which is then followed by the FSI transport (stage 4). Stage 4 conserves probability, so it does not affect the normalization.

## 2. Factorization into stages

Define the per-channel total cross section
```math
\sigma_{c,A}(E) = \int d\mathcal{I}_A\, d\mathcal{K}\, \frac{d\sigma_c}{d\mathcal{K}}
```
and the target-unit total $`\sigma_\text{tot}(E)=\sum_{A,c} n_A\,\sigma_{c,A}(E)`$. Eq. (1) then factorizes exactly:
```math
p = \underbrace{\frac{\Phi(E)\,\sigma_\text{tot}(E)}{\int \Phi\,\sigma_\text{tot}}}_{E_\nu\ \text{sampling}}
\;\times\;\underbrace{\frac{n_A\,\sigma_{c,A}(E)}{\sigma_\text{tot}(E)}}_{\text{stage 1}}
\;\times\;\underbrace{\frac{d\mathcal{I}_A\,\dfrac{d\sigma_c}{d\mathcal{K}}\,d\mathcal{K}}{\sigma_{c,A}(E)}}_{\text{stages 2 + 3}} .
```

**Stages 2 and 3.**
- Stage 2 (`NuclearModel::Sample`) draws $`\mathcal{I}`$ from the nuclear model alone.
- Stage 3 (`KinematicsGenerator::Generate`) draws $`\mathcal{K}`$ from a proposal density $`g(\mathcal{K}\mid E,\mathcal{I})`$ and returns the weight
```math
w = \frac{1}{g}\,\frac{d\sigma_c}{d\mathcal{K}}(E,\mathcal{I},\mathcal{K}) .
```
- The expectation of $`w`$ over stages 2 and 3 is the cross section per struck nucleon (or pair), $`\langle w\rangle = \sigma_{c,A}(E)/m_{c,A}`$. Here $`m_{c,A}`$ is the channel multiplicity, e.g. the number of neutrons for $`\nu`$ CCQE.
- The generator accepts a trial with probability $`w/w_\text{env}(E)`$ and otherwise repeats stages 2 and 3 **at fixed $`E`$, $`A`$, $`c`$**. The accepted trials therefore follow $`d\mathcal{I}_A\,\frac{d\sigma_c}{d\mathcal{K}}\,d\mathcal{K}/\sigma_{c,A}`$ exactly, provided the envelope is valid ($`w\le w_\text{env}`$).

**Stage 1 and energy sampling** need $`\sigma_{c,A}(E)`$, which `XSecTableBuilder` tabulates by Monte Carlo integration of the **same** weight:
```math
\sigma_{c,A}(E_i) = m_{c,A}\,\frac{1}{N}\sum_{j=1}^{N} w_j ,\qquad w_\text{max}(E_i) = \max_j w_j ,
```
evaluated at log-spaced nodes $`E_i`$. Between nodes, $`\sigma`$ is interpolated linearly. The envelope at energy $`E`$ is $`w_\text{env}(E) = f_\text{safety}\cdot\max\big(w_\text{max}(E_i), w_\text{max}(E_{i+1})\big)`$, where $`E_i \le E < E_{i+1}`$ and $`f_\text{safety}=1.2`$ by default. Because the tables and the events share one weight function, the totals and the differential distributions are consistent by construction. Adding a channel needs no channel-specific integration code.

**Energy sampling** (`EnergySampler`):
- For a histogram flux, $`\Phi`$ is constant inside each bin.
  - A bin $`i`$ is chosen with probability $`\propto \Phi_i\int_\text{bin}\sigma_\text{tot}\,dE`$. The integral is exact for the piecewise-linear $`\sigma_\text{tot}`$: trapezoids between the bin edges and all interior nodes.
  - $`E`$ is then drawn inside the bin by rejection against the bin maximum of $`\sigma_\text{tot}`$.
- For a mono-energetic beam, the table has a single node at the beam energy.

## 3. Envelope violations

If a trial has $`w > w_\text{env}`$, the event is kept with probability one (the least biased choice) and counted. The run summary reports two quantities:
- the number of such events;
- the largest $`w/w_\text{env}`$ seen in any trial. A value $`< 1`$ means the envelope was never violated.

For the CCQE model the run summary of every test run so far reports zero violations (Sec. 6 of [quasielastic.md](quasielastic.md)).

## 4. Normalization of the output

Events are unweighted. The flux-averaged total cross section per target unit,
```math
\langle\sigma\rangle_\Phi = \frac{\int \Phi(E)\,\sigma_\text{tot}(E)\,dE}{\int \Phi(E)\,dE},
```
is stored in the `run` tree as `flux_avg_xsec` (10⁻³⁸ cm²). It is also written as the NuHepMC attribute `NuHepMC.FluxAveragedTotalCrossSection` (convention G.C.2, with `TargetScale` "PerTargetMolecule" for mixtures or "PerAtom" for a single nucleus).

A flux-averaged differential cross section per target unit in a bin of width $`\Delta X`$ is
```math
\left\langle \frac{d\sigma}{dX}\right\rangle_\Phi = \langle\sigma\rangle_\Phi \;\frac{N_X}{N_\text{events}\,\Delta X}.
```

Two conversions are common:
- **Per nucleon:** divide by $`\sum_A n_A A`$ (13 for CH).
- **Per nucleus $`A`$ only:** select events with `target_pdg == A` and divide by $`n_A`$.

## 5. Reproducibility and threading

- **Random streams:** every batch of `batch_size` events uses its own xoshiro256++ stream, seeded through SplitMix64 from (seed, batch index). Every cross-section node is seeded the same way, from (table seed, node index).
- **Ordering:** batches are written strictly in order by a single writer thread.
- **Result:** the output is a deterministic function of the configuration and seed, independent of the number of threads (`tests/test_generator.cc`).
- **Thread safety:** no global mutable state is used in the hot path; models are immutable after construction.

## 6. Interfaces for future channels

- `HitNucleonRequest` already carries the number and species of nucleons a channel needs. A 2p2h/MEC channel requests two nucleons from the same `NuclearModel`. Its `Multiplicity` is the number of pairs, or 1 if the model's normalization already includes the pair counting.
- The event record separates the primary vertex (status 22 for hadrons that enter FSI) from the FSI output (status 1). An intranuclear cascade only needs to consume status-22 particles and append status-1 ones.
