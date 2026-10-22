# Second resonance region: N(1440), N(1520), N(1535) (CCNSTAR)

This note covers the excitation of the three isospin-1/2 resonances of the second resonance region by
the charged current,
$$
\nu_\ell\,n\to\ell^-R^+,\qquad \bar\nu_\ell\,p\to\ell^+R^0,\qquad R=P_{11}(1440),\ D_{13}(1520),\ S_{11}(1535),
$$
followed by $R\to N\pi$, $N\eta$ or $N\pi\pi$. It also covers the companion change in `CCDIS`: the
single-pion part of the inelastic continuum above $W=1.6$ GeV is rescaled per isospin channel by a fit to
deuterium data that includes these resonances (see [dis.md](dis.md) §3).

| Stage | Class | Role |
|---|---|---|
| 1 | `channel/CCNStarChannel` | $\nu$: neutrons only, $\bar\nu$: protons only (the $W^\pm$ cannot make an $I=1/2$ state from the other nucleon); multiplicity $N$ or $Z$; NuHepMC process id **401** (CC resonant) |
| 2 | `initial_state/RelativisticFermiGas`, `FreeNucleon` | struck nucleon (as CCRES) |
| 3 | `kinematics/NStarKinematics` | lepton kinematics, weight, resonance choice and decay |

**Why.** Below $W=1.6$ GeV, v4 had only the Δ(1232) plus the HNV chiral background (CCRES); above it, the
Bodek–Yang inclusive continuum (CCDIS). This had three consequences:
- On the neutron the $I=1/2$ second-resonance peak at $W\approx1.5$ GeV was missing.
- $d\sigma/dW$ stepped up by ×3.2 at 1.6 GeV.
- The duality-averaged continuum just above 1.6 GeV produced too many hard single-pion (especially $p\pi^0$) final states.

All three ingredients come from the theory paper by Lalakulich, Paschos and Piranishvili [1] (LPP):
- vector form factors fitted by LPP to the JLab/MAID electroproduction helicity amplitudes;
- axial couplings from PCAC and the measured $R\to N\pi$ widths.

Nothing in this channel is fitted to the MINERvA corpus.

---

## 1. Cross section

For a struck nucleon $p$ (on shell, de Forest prescription as in CCRES) and a resonance four-momentum
$p_R=k+p_\text{bound}-k'$ with $W^2=p_R^2$,
$$
d\sigma=\frac{\overline{|\mathcal M|^2}}{4E_\nu E_N}\,d\Phi_2\,\mathcal A_R(W^2)\,dW^2,\qquad
\overline{|\mathcal M|^2}=\frac{G_F^2V_{ud}^2}{4}L_{\mu\nu}H^{\mu\nu},\qquad
d\Phi_2=\frac{|\mathbf p^*_\ell|}{16\pi^2\sqrt s}\,d\Omega^*_\ell .
$$
$L_{\mu\nu}H^{\mu\nu}$ carries full spin sums: the factor 1/4 contains $(1/\sqrt2)^2$ from the four-fermion
coupling and 1/2 from the average over the nucleon spin. The resonances are summed **incoherently**,
$d\sigma=\sum_R d\sigma_R$. There is no interference with each other or with the Δ + background of
CCRES; this is the LPP approximation and the one used by Leitner *et al.* [3].

**Spectral function.** A relativistic Breit–Wigner with an energy-dependent width:
$$
\mathcal A_R(W^2)=\frac1\pi\,\frac{M_R\Gamma_R(W)}{(W^2-M_R^2)^2+M_R^2\Gamma_R(W)^2},\qquad
\Gamma_R(W)=\Gamma_R^0\Big[x_\pi\frac{\rho_{N\pi}(W)}{\rho_{N\pi}(M_R)}+x_\eta\frac{\rho_{N\eta}(W)}{\rho_{N\eta}(M_R)}+x_{\pi\pi}\frac{\rho_{N\pi\pi}(W)}{\rho_{N\pi\pi}(M_R)}\Big].
$$
- Two-body channels: $\rho(W)=p^{*2L+1}/W$, where $L$ is the orbital angular momentum of the $N\pi$ decay ($P_{11}$: 1, $D_{13}$: 2, $S_{11}$: 0) and $L=0$ for $N\eta$.
- $N\pi\pi$: $\rho(W)=W^{-1}\int dm\,p^*(W\to N+m)\,p^*(m\to\pi\pi)$, the three-body phase space.

### 1.1 Vertices (LPP Eqs. 4.2, 4.10, 4.15)

$q=p_R-p$, $Q^2=-q^2$, $\mu=m_N+M_R$, $i\sigma^{\nu\rho}q_\rho=-\tfrac12[\gamma^\nu,\slashed q]$.

**Spin 1/2:** $\langle R|J^\nu|N\rangle=\bar u(p_R)\Gamma^\nu u(p)$, with the vector and axial parts
$$
V^\nu=\frac{g_1^V}{\mu^2}(Q^2\gamma^\nu+\slashed q q^\nu)+\frac{g_2^V}{\mu}i\sigma^{\nu\rho}q_\rho,\qquad
A^\nu=g_1^A\gamma^\nu+\frac{g_3^A}{m_N}q^\nu .
$$
The parity of the resonance fixes where $\gamma_5$ sits:
- $P_{11}$: $\Gamma^\nu=V^\nu-A^\nu\gamma_5$;
- $S_{11}$: $\Gamma^\nu=V^\nu\gamma_5-A^\nu$.

**Spin 3/2 ($D_{13}$):** $\langle R|J^\nu|N\rangle=\bar\psi_\lambda(p_R)d^{\lambda\nu}u(p)$ with ($C_3^A=C_4^A=C_6^V=0$, as in LPP)
$$
d^{\lambda\nu}=g^{\lambda\nu}\Big[\frac{C_3^V}{m_N}\slashed q+\frac{C_4^V}{m_N^2}p_R\!\cdot\!q+\frac{C_5^V}{m_N^2}p\!\cdot\!q\Big]
-q^\lambda\Big[\frac{C_3^V}{m_N}\gamma^\nu+\frac{C_4^V}{m_N^2}p_R^\nu+\frac{C_5^V}{m_N^2}p^\nu\Big]
+\Big[g^{\lambda\nu}C_5^A+q^\lambda q^\nu\frac{C_6^A}{m_N^2}\Big]\gamma_5 .
$$
For positive parity (a $P_{33}$-like state) the code moves $\gamma_5$ to the vector part instead.

**Hadron tensors.** Spin 1/2:
$$
H^{\mu\nu}=\mathrm{Tr}[(\slashed p_R+M_R)\Gamma^\mu(\slashed p+m)\bar\Gamma^\nu] .
$$
Spin 3/2:
$$
H^{\mu\nu}=\mathrm{Tr}[P_{\rho\lambda}(p_R)\,d^{\lambda\mu}(\slashed p+m)\,\bar d^{\rho\nu}],
$$
with the Rarita–Schwinger spin sum
$$
P^{\rho\lambda}=-(\slashed p_R+M_R)\Big[g^{\rho\lambda}-\tfrac13\gamma^\rho\gamma^\lambda-\tfrac{2p_R^\rho p_R^\lambda}{3M_R^2}+\tfrac{p_R^\rho\gamma^\lambda-p_R^\lambda\gamma^\rho}{3M_R}\Big] .
$$
This is the projector used for the Δ in `PionProductionCurrent`.

Both tensors are evaluated numerically with the explicit Dirac matrices of `core/Dirac.h` and contracted with
the lepton tensor exactly as in CCRES. The closed-form structure functions of LPP serve as an
independent check (§4).

### 1.2 Form factors (LPP §IV, Table I)

$D_V=(1+Q^2/M_V^2)^2$ with $M_V=0.84$ GeV; $D_A=(1+Q^2/M_A^2)^2$ with $M_A=1.05$ GeV. Logarithms are $\ln(1+Q^2/1\,\text{GeV}^2)$.

**Vector form factors.** The weak (isovector) form factors of an $I=1/2$ state are $C^V=C^n-C^p$ (LPP Eq. 3.3).
For $P_{11}$ and $S_{11}$ the isoscalar part is neglected ($A^n=-A^p$), so $g^V_i=-2g_i^{(p)}$.

| | proton e.m. form factors (fit to helicity amplitudes) |
|---|---|
| $P_{11}(1440)$ | $g_1^{(p)}=\dfrac{2.3/D_V}{1+Q^2/4.3M_V^2}$, $g_2^{(p)}=-\dfrac{0.76}{D_V}[1-2.8\ln(1+Q^2)]$ |
| $S_{11}(1535)$ | $g_1^{(p)}=\dfrac{2.0/D_V}{1+Q^2/1.2M_V^2}[1+7.2\ln(1+Q^2)]$, $g_2^{(p)}=\dfrac{0.84}{D_V}[1+0.11\ln(1+Q^2)]$ |
| $D_{13}(1520)$ | $C_3^{(p)}=\dfrac{2.95/D_V}{1+Q^2/8.9M_V^2}$, $C_4^{(p)}=\dfrac{-1.05/D_V}{1+Q^2/8.9M_V^2}$, $C_5^{(p)}=-\dfrac{0.48}{D_V}$; $C_3^{(n)}=\dfrac{-1.13/D_V}{1+Q^2/8.9M_V^2}$, $C_4^{(n)}=\dfrac{0.46/D_V}{1+Q^2/8.9M_V^2}$, $C_5^{(n)}=-\dfrac{0.17}{D_V}$ |

**Axial form factors** (PCAC with the $R\to N\pi$ coupling, LPP App. A):
- $P_{11}$: $g_1^A(0)=-0.51$, $g_3^A=\dfrac{(M_R+m_N)m_N}{Q^2+m_\pi^2}g_1^A$;
- $S_{11}$: $g_1^A(0)=-0.21$, $g_3^A=\dfrac{(M_R-m_N)m_N}{Q^2+m_\pi^2}g_1^A$;
- $D_{13}$: $C_5^A(0)=-2.1$, $C_6^A=\dfrac{m_N^2C_5^A}{m_\pi^2+Q^2}$.

All of them use the shape $(\cdot)(0)/[D_A(1+Q^2/3M_A^2)]$. LPP chose the signs so that $W_3>0$ for every
resonance ($\sigma_\nu>\sigma_{\bar\nu}$). The test in §4 confirms that the numerical tensor reproduces their
$W_3$, sign included.

**Masses, widths and branching fractions.** $M_R$, $\Gamma^0_R$ and the elasticity $x_\pi$ are LPP Table I. The
axial couplings were derived from these values, so they are kept together. The inelastic remainder is taken
from the PDG 2024 ranges [4]:

| | $M_R$ | $\Gamma^0_R$ | $L_{N\pi}$ | $x_\pi$ | $x_\eta$ | $x_{\pi\pi}$ |
|---|---|---|---|---|---|---|
| $P_{11}(1440)$ | 1.440 | 0.350 | 1 | 0.60 | 0 | 0.40 |
| $D_{13}(1520)$ | 1.520 | 0.125 | 2 | 0.50 | 0 | 0.50 |
| $S_{11}(1535)$ | 1.535 | 0.150 | 0 | 0.40 | 0.45 | 0.15 |

### 1.3 W range and the boundary with CCDIS

The resonances are generated for $W<w_\text{max}=1.6$ GeV, the CCRES `w_max` and CCDIS `w_min`. Above 1.6 GeV
the Bodek–Yang structure functions describe the duality-averaged inclusive strength, which already contains
the resonance tails, so the line shapes are truncated rather than double counted. For a free nucleon the
fraction of $\int\mathcal A\,dW^2$ below 1.6 GeV is 0.46 ($P_{11}$), 0.62 ($D_{13}$) and 0.56 ($S_{11}$).

## 2. Decays

The resonance is chosen with probability $\propto d\sigma_R$ at the sampled point and the decay channel with
probability $\propto\Gamma_i(W)$. The charge states follow isospin (Clebsch–Gordan; quoted for $R^+$, while
$R^0$ is the mirror $p\leftrightarrow n$, $\pi^+\leftrightarrow\pi^-$):
- $N\pi$: $n\pi^+$ 2/3, $p\pi^0$ 1/3;
- $N\eta$: $p\eta$ (the η is left undecayed, status 1, like the π⁰);
- $N\pi\pi$ through $\pi\Delta$ in $I=1/2$, i.e. $|\tfrac12,\tfrac12\rangle=\sqrt{\tfrac12}\Delta^{++}\pi^-
  -\sqrt{\tfrac13}\Delta^+\pi^0+\sqrt{\tfrac16}\Delta^0\pi^+$ followed by $\Delta\to N\pi$:
  $p\pi^+\pi^-$ 5/9, $p\pi^0\pi^0$ 2/9, $n\pi^+\pi^0$ 2/9.

**Kinematics.** Decays are isotropic in the resonance rest frame, with two- or three-body phase space
(`DISHadronization::PhaseSpace`). The decay nucleon is Pauli blocked as in CCRES. The hadrons are status 22
(mother: the struck nucleon) and go through the cascade like CCRES hadrons, with no formation zone.

## 3. Sampling

- **W proposal:** 25% flat on $[M_N+m_{\pi^0},\min(w_\text{max},\sqrt s-m_\ell)]$; otherwise a truncated Cauchy around one of the resonances (equal shares).
- **Lepton angle:** the CCRES proposal, in the CM frame of $k+p_\text{bound}$.

**Cost on the ME flux.** The ¹²C table (150 nodes × 2×10⁵ samples) takes 55 s on 16 threads; generation is negligible.

## 4. Validation (`tests/test_nstar.cc`)

1. **Hadron tensor vs LPP closed forms.** For a nucleon at rest, a massless lepton, $W=M_R$, $E_\nu$ = 1.5, 3 and 8 GeV and $Q^2$ = 0.05, 0.4 and 1.2 GeV², both helicities, the numerical $L_{\mu\nu}H^{\mu\nu}$ equals
   $$32m_N\Big\{W_1Q^2+\tfrac{W_2}{m_N^2}\big[2(k\cdot p)(k'\cdot p)-\tfrac12m_N^2Q^2\big]\pm\tfrac{W_3}{m_N^2}\big[Q^2k\cdot p-\tfrac12p\cdot qQ^2\big]\Big\},$$
   with LPP's $W_i=(1/m_N)V_i$ (spin 1/2, Eqs. 4.11–4.13) or $(2/3m_N)V_i$ (spin 3/2, Eqs. 4.6–4.8, lower sign).
   - The agreement is to $10^{-8}$ for all three resonances.
   - This checks the vertex structure, the RS projector, the normalisation (the same $\overline{|\mathcal M|^2}$ convention as CCRES) and the sign of the $V$–$A$ interference.
2. **Widths:** $\Gamma(M_R)=\Gamma^0_R$, $\Gamma_{N\pi}(M_R)=x_\pi\Gamma^0_R$, closed channels vanish below threshold, and $\int_{W<2}\mathcal A\,dW^2\in[0.6,1]$.
3. **Events** (C and CH, ν and ν̄, 3 GeV):
   - exact four-momentum and charge conservation;
   - one lepton, one nucleon and ≥ 1 meson;
   - $W\le1.6$ GeV and process id 401;
   - ν hits only neutrons, ν̄ only protons;
   - both η and two-pion final states occur.

**Free-neutron cross section** ($W<1.6$, ×10⁻³⁸ cm²):

| $E_\nu$ [GeV] | 1.0 | 1.5 | 1.9 | 2.7 | 4.5 | 10 |
|---|---|---|---|---|---|---|
| total | 0.024 | 0.11 | 0.15 | 0.20 | 0.24 | 0.27 |
| $n\pi^+$ / $p\pi^0$ | 0.010/0.005 | 0.040/0.020 | 0.057/0.028 | 0.073/0.037 | 0.087/0.044 | |
| $p\pi^+\pi^-$ | 0.005 | | 0.031 | 0.038 | 0.044 | |
| $p\eta$ | 0.000 | | 0.014 | 0.022 | 0.029 | |

- **$N\pi$ line shape.** On the neutron at 2.7 GeV, $d\sigma(N\pi)/dW$ is 0.55–0.62 ×10⁻³⁸ cm²/GeV at $W=1.48$–1.53 GeV. That is the size of the second peak in LPP's single-pion $d\sigma/dW$ on the neutron (§V of [1]; about 0.8 including the Δ tail).
- **Two-pion data.** $\nu n\to\mu^-p\pi^+\pi^-$ is measured as 0.03–0.10 at $E_\nu$ = 1.1–1.7 GeV (BNL [6], ANL [7], as shown by Nakamura *et al.* [5]). The N* gives 0.02–0.03 here, against ≤ 0.02 before. The model is still low in this channel (no non-resonant ππN below 1.6 GeV; see §6).

**Single-pion data and the CCDIS rescaling.** The ANL/BNL deuterium $p\pi^+$, $n\pi^+$ and $p\pi^0$ cross sections without a W cut were used:
- $p\pi^+$: Wilkinson *et al.* [2], Fig. 6;
- $n\pi^+$, $p\pi^0$: the same reanalysis technique applied to those channels (Rodrigues *et al.* [8]). Only the extracted data points are used, digitized from the figures to about ±5%.

Those cross sections are fitted with CCRES + CCNSTAR + $R_{3/2}S_{3/2}+R_{1/2}S_{1/2}$. Here $S_I$ are the isospin components of the CCDIS single-pion strength ([dis.md](dis.md) §3).
- Fitted values: $R_{3/2}=0.28\pm0.15$ and $R_{1/2}=0.68\pm0.11$ (correlation −0.18), with $\chi^2=11.1/14$.
- At $R=1$ the $\chi^2$ is 49.5.
- Without the N* the fit gives $R_{1/2}=1.03\pm0.11$ ($\chi^2=13.8/14$). The N* therefore takes over a third of the $I=1/2$ single-pion strength, which CCDIS had placed above 1.6 GeV.

**Inclusive** (MINERvA ME flux on CH, per nucleon):
- CCNSTAR has $\langle\sigma\rangle=0.097\times10^{-38}$ cm² (2.5% of the total).
- $d\sigma/dW$ [×10⁻³⁸ cm²/GeV] is 0.91, 1.07, 1.67 and 1.69 in 1.4–1.5, 1.5–1.6, 1.6–1.7 and 1.7–1.8 GeV. The step at 1.6 GeV drops from ×3.2 to ×1.56.
- $\sigma_\text{CC}/E$ on carbon is 0.667 (4.5 GeV) and 0.661 (6.5 GeV), against MINOS 0.711 ± 0.033 and 0.708 ± 0.032 [9].

## 5. Configuration

```yaml
- type: CCNSTAR
  nuclear_model: {type: RFG}
  w_max: 1.6                     # = CCRES w_max = CCDIS w_min
  MV: 0.84                       # LPP vector mass
  MA: 1.05                       # LPP axial mass
  resonances: [P11_1440, D13_1520, S11_1535]   # default: all three
  proposal_q2_scale: 1.0
```

## 6. Limitations

- **Incoherent sum.** The N* do not interfere with each other, with the Δ, or with the HNV background. In particular the background of CCRES already contains some $I=1/2$ non-resonant $N\pi$ strength near 1.4–1.6 GeV, so a part of the N* strength may be double counted.
- **Axial couplings** come from PCAC with one form factor per resonance ($C_3^A=C_4^A=0$ for the $D_{13}$, as in LPP); their $Q^2$ dependence is a modified dipole with $M_A=1.05$ GeV, not measured.
- **Decays are isotropic,** and $N\pi\pi$ is pure $\pi\Delta$ with three-body phase space; the $N\sigma$ and $N\rho$ components and the Δ peak in the Dalitz plot are neglected.
- **No $I=3/2$ resonances** ($\Delta(1600)$, $\Delta(1620)$, $\Delta(1700)$). These are the only higher states reachable from the proton by $\nu$.
- **No non-resonant $\pi\pi N$ below 1.6 GeV.** ν n → p π⁺ π⁻ remains below the ANL/BNL data.
- **Truncation at $W=1.6$ GeV.** The tails above it are represented only by the duality-averaged Bodek–Yang strength.
- **No in-medium width or shift** of the resonances (as for the Δ).

## References

1. O. Lalakulich, E.A. Paschos, G. Piranishvili, "Resonance production by neutrinos: The second resonance region", *Phys. Rev. D* **74**, 014009 (2006) [hep-ph/0602210].
2. C. Wilkinson *et al.*, *Phys. Rev. D* **90**, 112017 (2014) [arXiv:1411.4482].
3. T. Leitner, O. Buss, L. Alvarez-Ruso, U. Mosel, *Phys. Rev. C* **79**, 034601 (2009) [arXiv:0812.0587].
4. S. Navas *et al.* (Particle Data Group), *Phys. Rev. D* **110**, 030001 (2024).
5. S.X. Nakamura, H. Kamano, T. Sato, *Phys. Rev. D* **92**, 074024 (2015) [arXiv:1506.03403], Fig. 19.
6. T. Kitagaki *et al.* (BNL), *Phys. Rev. D* **34**, 2554 (1986).
7. D. Day *et al.* (ANL), *Phys. Rev. D* **28**, 2714 (1983).
8. P. Rodrigues, C. Wilkinson, K. McFarland, *Eur. Phys. J. C* **76**, 474 (2016) [arXiv:1601.01888], App. A (data reanalysis only).
9. P. Adamson *et al.* (MINOS), *Phys. Rev. D* **81**, 072002 (2010) [arXiv:0910.2201].
