# 12C shell-model spectral function with short-range-correlated pairs (`SF`)

| Stage | Class | Role |
|---|---|---|
| 2 | `initial_state/ShellModelSpectralFunction` | Joint sampling of the struck-nucleon momentum $\mathbf k$, removal energy $E_R$, interaction point $\mathbf r$ and (for correlated pairs) the emitted partner nucleon |
| 3 | `kinematics/QuasiElasticKinematics` | Unchanged impulse-approximation vertex (de Forest), local Pauli blocking, partner written as a pre-FSI hadron |
| 4 | `fsi/IntranuclearCascade` | Starts all vertex hadrons from the stage-2 interaction point |

Configuration: `nuclear_model: {type: SF}` in a `CCQE` channel (12C only). It replaces the global
relativistic Fermi gas (RFG, [quasielastic.md](quasielastic.md) §4) of the primary CCQE vertex. The
other channels (CCRES, CCNSTAR, CCDIS, CCMEC) keep the RFG; `ModelFactory` refuses `SF`/`LFG` for them.

The construction follows the logic of Ankowski and Sobczyk [1] (mean-field shells plus a correlated
part) and of Benhar *et al.* [2,3] (impulse approximation with a realistic $P(\mathbf k,E)$), but every
ingredient was built here from the cited data; no generator code or tuned table was used.

---

## 1. Formalism

In the impulse approximation the CCQE cross section on 12C is the incoherent sum over nucleons of the
free-nucleon cross section weighted by the hole spectral function $P(\mathbf k, E)$, the probability
of removing a nucleon of momentum $\mathbf k$ leaving the residual system with excitation $E-S_N$:
$$
P(\mathbf k,E) = Z_\text{MF}\sum_{a=1s,1p}\frac{n_a}{N_a}\,|\phi_a(\mathbf k)|^2\,L_a(E)
\;+\;(1-Z_\text{MF})\,P_\text{SRC}(\mathbf k,E),\qquad \int P\,d^3k\,dE=1 .
$$

### 1.1 Mean field
- **Orbitals.** Harmonic-oscillator 1s and 1p orbitals with the oscillator length of the 12C point
  density used by the cascade, $b=1.590$ fm (fixed by the measured charge radius, [fsi_cascade.md](fsi_cascade.md)).
  In momentum space $|\phi_{1s}|^2\propto e^{-k^2b^2}$ and $|\phi_{1p}|^2\propto k^2b^2e^{-k^2b^2}$, i.e.
  $|\mathbf k|$ is a $\chi$ variable with 3 (1s) or 5 (1p) degrees of freedom and width
  $\sigma_k=\hbar c/(b\sqrt2)=87.7$ MeV per component. The interaction point is drawn from the same
  orbital in coordinate space ($\sigma_r=b/\sqrt2$).
- **Occupancies:** independent-particle shell model, 2 nucleons of each species in 1s and $Z-2=4$ in 1p.
- **Removal-energy distributions** $L_a(E)$: Lorentzians with the energy-dependent spreading width of
  Brown and Rho [4], $\Gamma(E)=a\,x^2/(b+x^2)$, $x=E-E_F$, $E_F=S_p=15.96$ MeV, $a=24$ MeV,
  $b=590$ MeV$^2$, truncated at $E_\text{max}=100$ MeV. Centroids from 12C$(e,e'p)$ [5,6]: 1s at 39 MeV;
  1p at 16.6 MeV (strength-weighted 3/2$^-$, 1/2$^-$, 3/2$^-$ states of 11B). The 1p strength is then
  almost a $\delta$ function at the Fermi level (width → 0 there), as observed. Neutron removal energies
  are shifted by $S_n-S_p=2.77$ MeV (Coulomb).
- **Mean removal energies** (sampled): protons 28.4 MeV, neutrons 31.2 MeV over the mean field (the
  Lorentzian tails skew the mean above the $\frac13E_{1s}+\frac23E_{1p}$ centroid average). The
  MINERvA STKI analysis reconstructs $p_n$ with $b=27.13$ MeV, close to the proton value.

### 1.2 Short-range-correlated pairs
- **Strength:** $1-Z_\text{MF}=0.20$ of the nucleons are in SRC pairs (12C$(e,e'pN)$, Subedi *et al.* [7]).
- **Momentum:** $k_\text{SRC}=0.25<k<k_\text{max}=1.0$ GeV with $n(k)\propto k^{-4}$ (contact/universal
  scaling; Hen *et al.* [8]), sampled exactly ($k^2n(k)\propto k^{-2}$).
- **Partner:** momentum $-\mathbf k+\mathbf K_\text{cm}$, $\mathbf K_\text{cm}$ Gaussian with
  $\sigma_\text{cm}=0.15$ GeV per component (Cohen *et al.* [9]); other isospin with probability 0.95
  (np dominance, Subedi [7], Hen *et al.* [10]).
- **Energy balance:** the partner is emitted on shell and the $A-2$ residual (10B for np, 10Be for pp,
  10C for nn) is left in its ground state with recoil $\mathbf K_\text{cm}$:
  $E_\text{hit}=M_A-E_2-\sqrt{M_{A-2}^2+K_\text{cm}^2}$, $M_{A-2}=M_A-M_1-M_2+S_{2N}$
  (two-nucleon thresholds from AME2020: $S_{np}=27.41$, $S_{pp}=27.19$, $S_{nn}=31.84$ MeV). The struck
  nucleon's removal energy $M-E_\text{hit}$ therefore grows with $k$ (≈ $k^2/2M$ + threshold), the
  $E$–$k$ correlation characteristic of the correlated part of realistic spectral functions [2,3].
- **Emission:** the partner is written as a status-22 (pre-FSI) hadron at the interaction point (drawn
  from the total density) and goes through the cascade. (`emit_partner: false` keeps it in the remnant
  instead.)

### 1.3 Energy convention, Pauli blocking and the vertex
- The struck nucleon is off shell with $E=M_A-\sqrt{M_{A-1}^{*2}+k^2}$, $M^*_{A-1}=M_A-M+E_R$ (residual
  recoil included); the outgoing nucleon is produced on shell with $\mathbf p'=\mathbf k+\mathbf q$ and
  energy $E+\omega$ (asymptotic convention, no outgoing potential at the vertex). The current is
  evaluated with the de Forest prescription ([quasielastic.md](quasielastic.md) §5.1).
- **Pauli blocking** uses the local Fermi momentum of the outgoing species at the interaction point,
  $k_F(r)=[3\pi^2\rho_{p,n}(r)]^{1/3}$ of the cascade density ($k_F(0)=0.273$ GeV).
- **Vertex → cascade:** the event stores the interaction point (`InitialState::position`) and the
  cascade starts every vertex hadron there, so the primary vertex and FSI see the same nucleus
  (previously the RFG vertex and the local-density cascade were independent; quasielastic.md §8).
- Four-momentum: beam + target = lepton + nucleon(s) + remnant exactly; the remnant carries $E_R$ as
  excitation (status 23), as before.

---

## 2. Implementation notes

- `Sample` draws shell / SRC, then $\mathbf k$, $\mathbf r$ and $E_R$; $L_a(E)$ is sampled from a
  20000-point inverse-CDF table built at construction.
- `PauliBlocked(…, InitialState)` and `Medium(…)` use `NuclearDensity::For(nucleus)` (thread-local
  cache shared with the cascade).
- `InitialState::hasSpectator/spectator` carry the SRC partner; `FillRemnant` removes it from the
  remnant; `QuasiElasticKinematics` writes it after the outgoing nucleon.
- Process id unchanged (CCQE, 200): the partner is part of the QE final state.
- No cache-version bump was needed: weights for the unchanged RFG YAML nodes are bit-identical.

## 3. Validation

All numbers from `build/nuflow-lfg-validate` (tables in `scratch/lfg/` of the step-17 workspace).

**Momentum distribution** ($k^2n(k)$, neutrons): $\langle k\rangle$ = 0.166 (RFG), 0.158 (LFG),
0.229 GeV (SF); 6.6% of SF nucleons have $k>0.5$ GeV (RFG/LFG: 0). The mean-field part reproduces
$\langle k^2\rangle_\text{HO}=(\frac13\cdot\frac32+\frac23\cdot\frac52)(\hbar c/b)^2$ to 1% (unit test). The
SF has no sharp edge at 0.221 GeV and, unlike the LFG of the same density, no excess at $k<0.08$ GeV
(the LFG piles up surface nucleons at small $k$, the known semiclassical artefact for light nuclei).

**Unit tests** (`tests/test_lfg.cc`, `[sf]`): SRC fraction and np share, residual nucleus species,
exact four-momentum balance of hit + partner + residual, shell medians at the centroids, mean removal
energy; end-to-end SF + cascade events conserve four-momentum to $10^{-9}$ GeV and start the cascade at
the stage-2 vertex.

**Inclusive 12C$(e,e')$ in the QE peak** (vector current, plane-wave IA, no FSI; data from the QES
archive, `data/electron/qes_12C.dat`), region $\omega\le\omega_\text{peak}+30$ MeV:

| $E$ [GeV], $\theta$ | data peak $\omega$ | RFG peak, height/data | LFG peak, height/data | SF peak, height/data |
|---|---|---|---|---|
| 0.480, 60° | 0.112 | 0.120, 1.07 | 0.100, 1.27 | 0.130, 0.89 |
| 0.560, 36° | 0.072 | 0.080, 1.23 | 0.060, 1.47 | 0.100, 0.94 |
| 0.680, 36° | 0.092 | 0.100, 1.12 | 0.080, 1.35 | 0.120, 0.90 |
| 0.961, 37.5° | 0.193 | 0.190, 1.05 | 0.170, 1.28 | 0.210, 0.90 |
| 1.299, 37.5° | 0.323 | 0.310, 0.97 | 0.290, 1.19 | 0.330, 0.84 |

The SF peak is broader and lower than the RFG peak (heights 0.84–0.94 of data vs 0.97–1.23), but in
PWIA it sits 10–30 MeV above the data, as is well known for spectral-function calculations without FSI:
the real part of the final-state optical potential shifts the peak down and the imaginary part
redistributes strength to the tails [2,3]. That correction is not included (§4). The RFG's position is
good because its $\epsilon_b=25$ MeV was historically fitted to these data.

**CCQE cross section** per neutron [$10^{-38}$ cm²], MA = 1.03 GeV (with the LQCD axial form factor used since step 25 the SF row becomes 0.929/1.023/0.994/0.980 at 1/3/6/10 GeV; [quasielastic.md](quasielastic.md) §7.2):

| $E_\nu$ [GeV] | 0.5 | 1 | 2 | 3 | 5 | 10 |
|---|---|---|---|---|---|---|
| RFG | 0.670 | 0.937 | 0.946 | 0.918 | 0.884 | 0.859 |
| SF | 0.494 | 0.764 | 0.795 | 0.777 | 0.754 | 0.732 |
| SF + RPA | 0.407 | 0.700 | 0.735 | 0.716 | 0.695 | 0.674 |

SF/RFG ≈ 0.85 above 1 GeV, roughly independent of $Q^2$. The mean field alone gives 0.945; the SRC
nucleons contribute about half the per-nucleon cross section of mean-field nucleons, because their
large off-shellness ($E_R\sim$ 0.1–0.5 GeV) moves the de Forest $\tilde Q^2$ above $Q^2$ and costs form
factor. The flux-averaged total of the MINERvA ME model falls from 50.26 to 49.35 ×10⁻³⁸ cm² per CH.

## 4. Limitations

- **No FSI on the inclusive response** (real optical potential / folding): the QE peak in energy
  transfer is 10–30 MeV too high in $(e,e')$ (§3). The cascade acts on hadrons only.
- **de Forest off-shell ambiguity** is large for SRC nucleons (factor ~2 per nucleon, §3); a
  current-conserving prescription would change the SRC share of the cross section.
- **Mean-field occupancies** are IPSM (1.0). The measured 1p spectroscopic strength (~0.6) is
  partly in the SRC part and partly at higher $E$ in the Lorentzian tails; the split is approximate.
- **Axial form factor:** the SF/RFG and SF/free ratios above were computed with the dipole. The SRC share
  depends on $F_A$ at the larger de Forest $\tilde Q^2$, so with the slower-falling LQCD $F_A$ (step 25)
  the SRC nucleons lose less: SF/free at 6 GeV rises from 0.76 to 0.78 (0.994/1.28).
- **12C only** (and ¹H via `FreeNucleon`). Parameters for other nuclei are not provided.
- **Long-range (RPA) correlations** are available (`rpa: true`, [local_fermi_gas.md](local_fermi_gas.md) §4), evaluated
  with the local density at the vertex. On top of the SF they are not used in the loop configuration
  (§5).
- CCMEC remains defined in the global RFG (Ruiz Simo *et al.*), so the 1p1h and 2p2h ground states
  differ; CCRES/CCNSTAR/CCDIS also keep the RFG.

## 5. RPA on top of the SF (Sec. 5.3 of the code comments)

`rpa: true` multiplies the 1p1h hadron tensor by the Valencia RPA factors of
[local_fermi_gas.md](local_fermi_gas.md) §4 at the local density $\rho(r)$ of the vertex. Two regularisations were
needed to keep the factors finite for spectral-function kinematics (both are documented parameters):
1. **Energy argument.** All of $U_N$, $U_\Delta$, $V_l$, $V_t$ are evaluated at the particle-hole energy
   $\nu=(|\mathbf k+\mathbf q|^2-k^2)/2M$, clamped into the open 1p1h support of the local Fermi sea,
   $[\max(0,q^2-2qk_F),\,q^2+2qk_F]/2M$ (a mean-field nucleon with $k>k_F(r)$ would otherwise give a
   real $U_N$ and spurious zeros of $1-UV$; evaluating $V_l$ at the large de Forest $\tilde q^0$ of a
   deeply bound nucleon reaches the pion pole).
2. **Widths.** The forward $\Delta h$ propagator carries the free Δ width, $1/(q^0-\omega_\Delta+i\Gamma_\Delta/2)$,
   $\Gamma_\Delta=117$ MeV (PDG), since the static zero-width bubble has a pole at
   $q^0=\sqrt{M_\Delta^2+q^2}-M$ that the 1p1h support reaches above $|\mathbf q|\approx0.6$ GeV; and $U_N$
   is evaluated at $q^0+i\Gamma_{ph}/2$ with $\Gamma_{ph}=5$ MeV, of the order of the spreading widths of
   the giant resonances of light nuclei, which damps the zero-sound pole of undamped nuclear matter
   just above the continuum at $|\mathbf q|<0.1$ GeV. With these the factors stay below 2.7 everywhere
   (they reached $10^{10}$ before); at $|\mathbf q|>0.15$ GeV the change from $\Gamma_{ph}$ is < 3%.

RPA reduces $\sigma_\text{CCQE}$ by another 8% (table above) and $d\sigma/dQ^2$ by 25–30% at
$Q^2<0.1$ GeV², crossing 1 at $Q^2\approx0.4$ GeV². In the step-17 corpus comparison (SF+RPA trial vs
the SF model v6) it flattened the QE-like $p_T$ shape below 0.85 GeV and improved ccqelike_3d_le
(χ² 4609 → 3738) and the STKI $\delta\alpha_T$/$\delta\phi_T$ shapes, but pushed ccinc_eavail_q3 at
$q_3=0.2$–0.6 GeV and $E_\text{avail}<0.06$ GeV down to 0.4–0.6 of the data (per-$q_3$ integrals
0.86/0.84/0.82 → 0.62/0.61/0.67 for $q_3<0.4$; χ² 103k → 157k) and the QE-like normalisation from 0.80
to 0.75. Since the SF's local Pauli blocking and removal-energy spread already provide the
low-$q_3$ reduction that ccinc_eavail_q3 asks for, and the RPA of Ref. [11] is derived for the LFG
particle-hole response rather than for a spectral function, the loop configuration leaves it off.

## 6. Configuration reference

```yaml
channels:
  - type: CCQE
    nuclear_model:
      type: SF            # 12C shell-model spectral function (this note); LFG: local_fermi_gas.md
      mf_fraction: 0.80   # Z_MF
      e1s: 0.039          # 1s centroid [GeV] (proton removal)
      e1p: 0.0166         # 1p centroid [GeV]
      width_a: 0.024      # Brown-Rho width a [GeV]
      width_b: 590e-6     # Brown-Rho width b [GeV^2]
      e_max: 0.100        # upper end of L_a(E) [GeV]
      k_src: 0.25         # SRC tail [k_src, k_max] [GeV]
      k_max: 1.0
      sigma_cm: 0.15      # pair c.m. width per component [GeV]
      np_fraction: 0.95
      emit_partner: true
    form_factors: {vector: BBBA05, axial: dipole, MA: 1.03}   # loop config since step 25: LQCD z expansion,
                          # quasielastic.md §1.3
    rpa: false            # true, or a map (local_fermi_gas.md §7)
```

## References

1. A. M. Ankowski and J. T. Sobczyk, "Construction of spectral functions for medium-mass nuclei", *Phys. Rev. C* **77**, 044311 (2008), [arXiv:0711.2031](https://arxiv.org/abs/0711.2031).
2. O. Benhar, N. Farina, H. Nakamura, M. Sakuda and R. Seki, *Phys. Rev. D* **72**, 053005 (2005), [hep-ph/0506116](https://arxiv.org/abs/hep-ph/0506116).
3. O. Benhar, A. Fabrocini, S. Fantoni and I. Sick, *Nucl. Phys. A* **579**, 493 (1994).
4. G. E. Brown and M. Rho, *Nucl. Phys. A* **372**, 397 (1981) (energy dependence of the hole spreading width).
5. L. Lapikás, G. van der Steenhoven, L. Frankfurt, M. Strikman and M. Zhalov, "Transparency of 12C for protons", *Phys. Rev. C* **61**, 064325 (2000).
6. J. Mougey *et al.*, *Nucl. Phys. A* **262**, 461 (1976) (12C$(e,e'p)$ 1s and 1p removal energies).
7. R. Subedi *et al.*, "Probing cold dense nuclear matter", *Science* **320**, 1476 (2008).
8. O. Hen *et al.*, "Evidence for the contact-term / $1/k^4$ scaling of high-momentum tails", *Phys. Rev. C* **92**, 045205 (2015).
9. E. O. Cohen *et al.* (CLAS), "Center of mass motion of short-range correlated nucleon pairs", *Phys. Rev. Lett.* **121**, 092501 (2018).
10. O. Hen *et al.* (CLAS), "Momentum sharing in imbalanced Fermi systems", *Science* **346**, 614 (2014).
11. J. Nieves, J. E. Amaro and M. Valverde, *Phys. Rev. C* **70**, 055503 (2004), [nucl-th/0408005](https://arxiv.org/abs/nucl-th/0408005).
