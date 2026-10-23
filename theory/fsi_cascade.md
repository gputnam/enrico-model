# Stage 4: `IntranuclearCascade` (semi-classical INC for nucleons and pions)

`fsi/IntranuclearCascade` (config `fsi: {type: cascade}`) transports the hadrons from the primary vertex (status 22) through the nucleus.
- Hadrons that escape are appended with status 1, with `mother` pointing to their primary ancestor.
- The excited residual nucleus is appended with status 25 (`kStatusRemnantPostFSI`). It takes the four-momentum that balances the event.

Energy, momentum, charge and baryon number are therefore conserved exactly, event by event. Hydrogen targets are transparent, as before. All ingredients come from hadron–nucleon or hadron–nucleus physics, and **nothing is tuned to the neutrino corpus**.

Files:
- `fsi/NuclearDensity.*`: density, local Fermi momentum, vertex sampling.
- `fsi/HadronNucleonXSec.*`: free NN and πN cross sections, Oset–Salcedo Δ self-energy.
- `fsi/IntranuclearCascade.*`: the cascade.
- `apps/nuflow-fsi-validate`: hadron–nucleus and event-level validation.
- `tests/test_cascade.cc`: unit tests.

## 1. Nucleus

**Density.** For A ≤ 16 the density is the harmonic-oscillator shell-model point-nucleon density
$$\rho(r)=\rho_0\,(1+c\,x^2)\,e^{-x^2},\qquad x=r/b,\quad c=(Z-2)/3 ,$$
with protons and neutrons scaled by Z/A and N/A.
- The length b is fixed by the measured rms charge radius ($r_\text{ch}(^{12}\text{C})=2.4702$ fm, Angeli & Marinova 2013), after the nucleon finite size is removed:
$$\langle r^2\rangle_\text{pt}=r_\text{ch}^2-r_p^2-\tfrac NZ r_n^2-\tfrac{3}{4M^2},$$
with $r_p=0.8409$ fm and $r_n^2=-0.1161$ fm².
- For HO, $\langle r^2\rangle=b^2(6+15c)/(4+6c)$.
- For ¹²C this gives b = 1.590 fm, point rms 2.340 fm, $\rho(0)=0.179$ fm⁻³ and $k_F(0)=0.273$ GeV.
- The cascade region ends at $R_\text{max}=5.52$ fm, where $\rho=10^{-4}\rho(0)$.

Heavier nuclei use a two-parameter Fermi density with a = 0.52 fm and the same rms constraint. This is untested.

**Local Fermi gas.** The local Fermi momentum is $k_{F,N}(r)=(3\pi^2\rho_N(r))^{1/3}$. It is used for:
- Pauli blocking of every outgoing nucleon;
- the momenta of the Fermi-sea collision partners (uniform in the local Fermi sphere).

**Vertex.** The primary vertex is taken from the stage-2 nuclear model when it provides one (`InitialState::hasPosition`: the CCQE spectral function and local Fermi gas, [spectral_function.md](spectral_function.md)); otherwise it is drawn from $r^2\rho(r)$. All primary hadrons start there, including both 2p2h nucleons. The position is stored in the new `Particle::x` field (fm): the vertex for status 22, the exit point for status 1.

**Nucleon potential.** Nucleons move in a local well $V_N(r)=T_F(r)+S_N$, with $T_F=\sqrt{k_F^2+M^2}-M$ and the measured separation energies $S_p=15.96$ and $S_n=18.72$ MeV (AME2020).
- The conserved quantity along a trajectory is the asymptotic kinetic energy $E^*=T_\text{in}(r)-V_N(r)$.
- Collision kinematics use the in-medium momentum.
- A Fermi-sea nucleon knocked out at r therefore pays $V_N(r)-t\ge S_N$.
- A nucleon whose $E^*\le0$ after a collision is captured by the remnant.

**RFG / local-density inconsistency.** For CCQE with `SF`/`LFG` (the loop default since step 17) the vertex, local Pauli blocking and the cascade share the density. The other channels' primary vertex still uses the global RFG ($k_F=0.221$ GeV, $\epsilon_b=25$ MeV, [quasielastic.md](quasielastic.md)), and the primary nucleon's on-shell kinetic energy is taken as its $E^*$. Without collisions the cascade therefore reproduces the transparent result exactly.
- Mean removal energy of the well: $\langle V-t\rangle=\tfrac25T_F+S\approx 26$ MeV at $k_F=0.221$. This is close to $\epsilon_b$, but not identical.
- Pauli blocking at the primary vertex uses the global $k_F$; the cascade uses the local one.

## 2. Nucleons

**Collision rate.** For each target species t (p, n), a partner is drawn from the local Fermi sea. The rate per unit path is
$$\lambda^{-1}=\sum_t \rho_t(r)\,f_t\,\sigma^{tot}_{Nt}(p_\text{lab})\,v_\text{Møller}/v_1 .$$
- $f_t$ is the fraction of species-t nucleons still in the remnant (the density shape is kept).
- $p_\text{lab}$ is the projectile momentum in the partner's rest frame.
- A collision occurs in a step of 0.2 fm with probability $1-e^{-\Delta s/\lambda}$.
- If an outgoing nucleon lies inside the local Fermi sphere, the collision is Pauli blocked and the particle continues (Metropolis et al. 1958).
- A knock-out must leave at least one nucleon in the remnant. The remnant carries the binding energy; this only matters for full break-up of ¹²C in high-multiplicity DIS events.

**Cross sections.** Cugnon, L'Hôte, Vandermeulen, NIM B 111, 215 (1996), fits to the pp/np data, p in GeV/c. The elastic forms above 0.44 GeV/c are as quoted in Buss et al., Phys. Rep. 512 (2012) App. B.

| | p < 0.44 | 0.44–0.8 | 0.8–2 | > 2 |
|---|---|---|---|---|
| σ_el^pp | 34 (p/0.4)^-2.104 | 23.5+1000(p−0.7)^4 | 1250/(p+50) − 4(p−1.3)² | 77/(p+1.5) |
| σ_el^np | 6.3555 p^-3.2481 e^{−0.377 ln²p} | 33+196\|p−0.95\|^2.5 | 31/√p | 77/(p+1.5) |

Total cross sections:
- **pp:** the elastic form below 0.8 GeV/c; then $23.5+24.6/(1+e^{-(p-1.2)/0.1})$ up to 1.5; then $41+60(p-0.9)e^{-1.2p}$ up to 5; then $48+0.522\ln^2p-4.51\ln p$.
- **np:** the elastic form below 1.0 GeV/c; then $24.2+8.9p$ up to 2; then 42 mb.
- **nn:** taken equal to pp.
- Below p = 0.1 GeV/c the values are frozen.

**Elastic scattering.** The angular distribution is $d\sigma/dt\propto e^{Bt}$, with the Cugnon slopes:
- pp: $B=5.5p^8/(7.7+p^8)$ for p < 2 GeV/c, and $5.334+0.67(p-2)$ above;
- np: 0 below 0.225 GeV/c; $16.53(p-0.225)$ up to 0.6; $7.16-1.63p$ up to 1.6; then as pp.

The distribution is symmetrised forward/backward with weight 1/2, which approximates the np charge-exchange backward peak.

**Inelastic scattering** ($\sigma^{tot}-\sigma^{el}$) is modelled as NN → NΔ → NNπ:
- Isospin factors: pp → nΔ⁺⁺ : pΔ⁺ = 3:1; np → pΔ⁰ : nΔ⁺ = 1:1; nn is the mirror of pp.
- The Δ mass is drawn from a relativistic Breit–Wigner with the CCRES width, times two-body phase space.
- The production angle uses the same $e^{Bt}$ form. The Δ decays isotropically with Clebsch–Gordan charge branchings.
- Above about 2 GeV/c, multi-pion production is lumped into this single-pion channel.

## 3. Pions: Δ-hole model with the Oset–Salcedo Δ self-energy

A pion with lab kinetic energy $T_\pi$ meets a Fermi-sea nucleon at invariant mass W. The Δ propagator in the medium is (Oset & Salcedo, NPA 468, 631 (1987); Salcedo et al., NPA 484, 557 (1988))
$$G_\Delta=\big[W-M_\Delta+\tfrac i2\tilde\Gamma-\Sigma_\Delta\big]^{-1},\qquad -\mathrm{Im}\,\Sigma_\Delta=C_Q\Big(\tfrac{\rho}{\rho_0}\Big)^{\alpha}+C_{A2}\Big(\tfrac{\rho}{\rho_0}\Big)^{\beta}+C_{A3}\Big(\tfrac{\rho}{\rho_0}\Big)^{2\beta}.$$

**Self-energy parameters.**
- $\rho_0=0.17$ fm⁻³.
- Each coefficient is a quadratic $a x^2+bx+c$ in MeV, with $x=T_\pi/m_\pi$:

| | a | b | c |
|---|---|---|---|
| C_Q | −5.19 | 15.35 | 2.06 |
| C_A2 | 1.06 | −6.64 | 22.66 |
| C_A3 | −13.46 | 46.17 | −20.34 |
| α | 0.382 | −1.322 | 1.466 |
| β | −0.038 | 0.204 | 0.613 |

- The parametrisation is valid for 85–315 MeV. Outside that range it is frozen at the end points.
- The real part $\mathrm{Re}\,\Sigma_\Delta$ is neglected.

**Free width and Pauli blocking.**
- The free width is the CCRES/HNV P-wave width, $\Gamma(W)=\frac{1}{6\pi}(f^*/m_\pi)^2\frac{M}{W}q_\text{cm}^3$, with f* = 2.14 (0.115 GeV at the pole).
- $\tilde\Gamma=\Gamma\cdot F_\text{Pauli}$, where $F_\text{Pauli}$ is the fraction of isotropic Δ decays, boosted to the lab, whose nucleon lies above the local $k_F$. It is computed analytically per collision.

**Channel cross sections.** Formation uses the free coupling. The decay or collision branch c then has partial width $\Gamma_c$:
$$\sigma_c=I\,\frac{8\pi}{q_\text{cm}^2}\,\frac{(\Gamma/2)(\Gamma_c/2)}{(W-M_\Delta)^2+\Gamma_\text{tot}^2/4},\qquad \Gamma_\text{tot}=\tilde\Gamma+2(C_Q+C_{A2}+C_{A3}),$$
where:
- $\Gamma_c/2=\tilde\Gamma/2+C_Q$ for quasi-elastic scattering;
- $\Gamma_c/2=C_{A2}$ for two-body absorption and $C_{A3}$ for three-body absorption;
- $I$ is the πN→Δ isospin factor: 1 for π⁺p and π⁻n, 2/3 for π⁰N, 1/3 for π⁺n and π⁻p.

In the free limit this is the P₃₃ Breit–Wigner, with σ(π⁺p) = 190 mb at the peak. The rate per unit length is $\sum_t\rho_t f_t\sigma_t v_\text{Møller}/v_\pi$, which is the Oset–Salcedo probability $-\mathrm{Im}\,\Pi/q$.

**Above the Δ region.** The Δ term is faded out linearly between W = 1.35 and 1.45 GeV. It is replaced by an isospin-averaged non-resonant πN cross section of 30 mb (the PDG πN total cross sections are 25–45 mb in this region). This term is treated as elastic-like scattering with $B=7$ GeV⁻².

**Final states.**
- **Quasi-elastic scattering and charge exchange.** The Δ decays with Clebsch–Gordan branchings, which gives charge exchange, for example π⁺n → Δ⁺ → π⁰p (2/3). The angular distribution in the CM frame is $1+3\cos^2\theta$. Because $\tilde\Gamma$ is already in the rate, the decay is re-drawn until the nucleon is unblocked.
- **Two-body absorption** (πNN → NN):
  - The second nucleon is drawn with per-pair weight 1 for np pairs and 0.05 for like (T = 1) pairs. The isoscalar quasi-deuteron channel dominates by about 20 (³He data; see the review arXiv:2210.06114 and Ashery & Schiffer, ARNPS 36, 207 (1986)).
  - Only charge-allowed pairs are used. A π⁺ is therefore absorbed mainly on np, giving pp; a π⁰ on np gives np.
  - The available four-momentum (pion plus the bound nucleons) goes into 2 nucleons, isotropic in their CM frame. Pauli blocking is handled by re-drawing, up to 20 times.
  - The nucleons then cascade and pay their removal energy.
- **Three-body absorption.** As two-body absorption, plus a third nucleon drawn from the remnant composition, with three-body phase space.

**Pions** feel no potential, so for a pion $E^*=T$.

## 4. Formation zone

Hadrons from CCDIS (process 600) travel freely for $l=|\vec p|/\mu^2$ before they can interact. This is the SKAT parametrisation with $\mu^2=0.08\pm0.04$ GeV² (Baranov et al., SKAT, IHEP preprint 84-4 (1984); discussed by Golan et al., PRC 86, 015505 (2012)). A 1 GeV/c pion therefore travels 2.5 fm. CCQE, CCMEC and CCRES hadrons have no formation zone; the Δ lifetime is neglected.

## 5. Validation (`apps/nuflow-fsi-validate`, `tests/test_cascade.cc`)

**Free inputs** (`nuflow-fsi-validate xsec`):
- NN: σ_pp = 146 mb and σ_np = 446 mb at 0.2 GeV/c (T = 21 MeV).
- σ(π⁺p) = 16, 76, 190 and 114 mb at T_π = 85, 125, 165 and 245 MeV.
- At ρ₀ and 165 MeV: $C_Q,C_{A2},C_{A3}$ = 13.0, 16.3, 15.4 MeV, i.e. a collisional width of 89 MeV.

**Local mean free paths** (π⁺, 165 MeV, ¹²C centre): λ_abs = 1.7 fm and λ_QE = 1.3 fm. Golan et al. Fig. 2 shows similar values for the Oset model: about 0.6 fm⁻¹ for each process at the centre of Fe.

**π⁺–¹²C** (200k pions per energy on a disk of radius 5.52 fm; cascade "reaction" = at least one interaction), compared with Ashery et al., PRC 23, 2173 (1981), Table I (σ_reac = σ_tot − σ_el):

| T_π [MeV] | reac (data) | abs (data) | cx (data) | inel (data) |
|---|---|---|---|---|
| 85 | 160 (287) | 83 (109±20) | 15 (35±12) | 62 (143±26) |
| 125 | 278 (417) | 100 (166±26) | 33 (38±12) | 144 (213±33) |
| 165 | 340 (447) | 94 (194±36) | 48 (46±23) | 196 (207±33) |
| 205 | 326 (412) | 90 (157±37) | 48 (45±23) | 186 (210±51) |
| 245 | 293 (366) | 84 (95±32) | 43 (47±23) | 163 (224±30) |
| 315 | 245 (309) | 68 (64±27) | 38 (45±22) | 136 (200±22) |

- Charge exchange agrees with the data. Absorption agrees at 245–315 MeV, and is about 15% low at 85 MeV.
- **At the Δ peak (125–205 MeV) absorption is low by a factor of about 2.** The published NuWro implementation of the same Oset model (Golan et al. Fig. 4) peaks at about 140 mb, also below Ashery.
- The reaction cross section is 20–40% low. A cascade with *free* πN cross sections gives the same σ_reac (344 mb at 165 MeV), so this deficit is a limit of the semi-classical straight-line geometry for beam pions (wavelength/diffraction effects at the surface, λ̄ = 0.73 fm). It is not a local-rate effect, and it matters less for pions produced inside the nucleus.
- DUET (arXiv:1506.07783) measures σ_abs+cx = 198–225 mb at 106–187 MeV; the model gives 134–142 mb.
- π⁻ gives the same values within 2%.
- **No factor was fitted to these data.**

**Proton transparency** (proton started at a density-weighted point, isotropic; T = probability of no unblocked NN collision):
- T = 0.61 at T_p = 0.2 GeV, 0.56 at 0.5, 0.505 at 1.0–1.5, and 0.51–0.52 at 2–4 GeV.
- Data: Garrow et al., PRC 66, 044613 (2002) give $T_C$ = 0.548, 0.570 and 0.573 at Q² = 3.3, 6.1 and 8.1 GeV² (T_p ≈ 1.8, 3.3 and 4.3 GeV). Earlier JLab/NE18 values are about 0.6 at lower Q².
- The model is about 0.05 (10%) more absorptive than the data. This is the expected size of the neglected NN-correlation (hole) and in-medium cross-section reductions.

**Neutrino events** (`nuflow-fsi-validate events config/minerva_me_fhc_numu_CH_loop.yaml 1000000`):
- Exact conservation: |ΔE|, |Δp| < 4×10⁻¹³ GeV, ΔQ = ΔB = 0.
- CCRES events on carbon with a primary π⁺: absorbed 19.8%, charge-exchanged to π⁰ 5.2%, single π⁺ 69.6%, other 5.4%.
- Mean final-state proton and neutron multiplicities per channel:

| channel | ⟨n_p⟩ | ⟨n_n⟩ |
|---|---|---|
| CCQE | 1.18 | 0.31 |
| 2p2h | 2.10 | 0.75 |
| CCRES | 1.62 | 0.67 |
| CCDIS | 1.21 | 0.74 |

- An FSI pion (NN → NΔ) is produced in 5.8% of CCQE and 8.8% of 2p2h events.
- Generation time grows by about 30% (2M loop events in about 1 min of generation).

## 6. Limitations

- Δ-peak pion absorption is about 2× below Ashery, and σ_reac is low for beam pions (see above). There is no s-wave πN term (b₀, b₁, Im B₀) and no Re Σ_Δ. The Oset coefficients are frozen below 85 MeV.
- Free NN cross sections are used, with no in-medium reduction or NN correlations. The np angular distribution is symmetrised. NN inelasticity is single-Δ only.
- Trajectories are straight lines with no refraction, and protons have no Coulomb barrier. The potential only enters through $E^*$; its momentum is absorbed by the remnant.
- The RFG primary vertex of the non-QE channels and the local-density cascade are inconsistent (Sec. 1); CCQE with `SF`/`LFG` is consistent.
- There is no de-excitation of the remnant: no evaporation nucleons and no γ's. The status-25 remnant keeps the excitation energy.
- The formation zone applies to DIS only. Δ propagation in CCRES is neglected.
- Only the ¹²C density (HO) has been validated.

## 7. Configuration

```yaml
fsi: {type: cascade, step: 0.2, nucleon_rescattering: true, nn_pion_production: true,
      pion_quasielastic: true, pion_absorption: true, like_pair_absorption: 0.05,
      formation_zone: true, formation_mu2: 0.08, formation_processes: [600]}
```

## References

- J. Cugnon, D. L'Hôte, J. Vandermeulen, Nucl. Instrum. Meth. B 111, 215 (1996).
- O. Buss et al., Phys. Rept. 512, 1 (2012) [arXiv:1106.1344], App. B.
- E. Oset, L.L. Salcedo, Nucl. Phys. A 468, 631 (1987).
- L.L. Salcedo, E. Oset, M.J. Vicente-Vacas, C. Garcia-Recio, Nucl. Phys. A 484, 557 (1988).
- T. Golan, C. Juszczak, J.T. Sobczyk, Phys. Rev. C 86, 015505 (2012) [arXiv:1202.4197].
- D.S. Baranov et al. (SKAT), IHEP preprint 84-4 (1984).
- D. Ashery et al., Phys. Rev. C 23, 2173 (1981).
- K. Ieki et al. (DUET), Phys. Rev. C 92, 035205 (2015) [arXiv:1506.07783].
- K. Garrow et al., Phys. Rev. C 66, 044613 (2002) [hep-ex/0109027].
- I. Angeli, K.P. Marinova, At. Data Nucl. Data Tables 99, 69 (2013).
- AME2020: W.J. Huang et al., Chin. Phys. C 45, 030002 (2021).
- D. Ashery, J.P. Schiffer, Ann. Rev. Nucl. Part. Sci. 36, 207 (1986).
- N. Metropolis et al., Phys. Rev. 110, 185 (1958).
