# Charged-current coherent pion production (`CCCOH`, process 100)

Stage 1 `CCCOHChannel`, stage 2 `GroundStateNucleus`, stage 3 `CoherentPionKinematics` (+ `PionNucleusElastic`),
tests in `tests/test_coherent.cc`.

$$\nu_\mu + A \to \mu^- + \pi^+ + A_\text{g.s.},\qquad \bar\nu_\mu + A \to \mu^+ + \pi^- + A_\text{g.s.}$$

The neutrino scatters coherently on the whole nucleus, which stays in its ground state. No quantum numbers
are exchanged with the nucleus and the four-momentum transfer to it, $|t| = -(q-p_\pi)^2$, is small
($|t|\lesssim 1/b \sim 0.02$ GeV²). The final state is the lepton, one forward pion and the recoiling
nucleus, with no nucleons. The channel applies only to nuclei ($A\ge 4$). On a free proton the same final
state is part of resonant/non-resonant single-pion production (`CCRES`).

## 1. Model: PCAC with measured pion–nucleus elastic scattering (Berger–Sehgal 2009)

Adler's theorem (S.L. Adler, Phys. Rev. 135, B963 (1964)) relates forward neutrino scattering
($Q^2\to0$) to the pion–hadron cross section, because only the divergence of the axial current, i.e. the
pion field (PCAC), survives there. Kopeliovich and Marage (Int. J. Mod. Phys. A 8, 1513 (1993)) extended
it to $Q^2>0$ and kept the charged-lepton mass. Ch. Berger and L.M. Sehgal, *PCAC and coherent pion
production by low energy neutrinos*, Phys. Rev. D 79, 053003 (2009) [arXiv:0812.2653] (BS09) apply it to the
coherent process. BS09 Eq. 6:

$$\frac{d\sigma^{CC}}{dQ^2\,dy\,dt} = \frac{G_F^2\cos^2\theta_C f_\pi^2}{2\pi^2}\,\frac{E}{|\mathbf q|}\,uv
\left[\Big(G_A-\frac{1}{2}\frac{Q^2_\text{min}}{Q^2+m_\pi^2}\Big)^2+\frac{y}{4}\,(Q^2-Q^2_\text{min})\,\frac{Q^2_\text{min}}{(Q^2+m_\pi^2)^2}\right]
\frac{d\sigma(\pi A\to\pi A)}{dt}$$

- **Kinematic variables.** $y=\nu/E$, $\nu = E-E'$, $|\mathbf q|=\sqrt{\nu^2+Q^2}$, and $u,v=(E+E'\pm|\mathbf q|)/2E$.
  The factor $E\,uv/|\mathbf q|$ is the exact kinematic factor of the forward-lepton theorem. Rein and
  Sehgal (RS83) replaced it by its $Q^2=0$ value $(1-y)/y$. The exact factor lowers the cross section at
  low $E$ (BS09 Fig. 1a).
- **Axial form factor.** $G_A=m_A^2/(m_A^2+Q^2)$.
- **Lepton-mass (pion-pole) correction.** $Q^2_\text{min}=m_\mu^2\,y/(1-y)$.
  - The first term is helicity non-flip. The pion pole of the hadronic axial current interferes
    destructively with it at $Q^2\lesssim m_\pi^2$: this is the "Adler screening" of D. Rein and L.M.
    Sehgal, Phys. Lett. B 657, 207 (2007) [hep-ph/0606185].
  - The second term is helicity flip and vanishes in the forward direction.
  - With the correction the CC/NC ratio is 1.37 at 1 GeV, 1.74 at 4 GeV and 1.84 at 25 GeV (model values;
    BS09 Fig. 1b), against the naive $2\cos^2\theta_C$.
- **Transverse contributions.** $\sigma_{L}$ and $\sigma_{R}$ vanish at $Q^2=0$ and are neglected, as in BS09.
  BS09 estimate the total theory error at about 20%.
- **Elastic cross section.** $d\sigma(\pi^+A\to\pi^+A)/dt$ is evaluated at the pion kinetic energy
  $T_\pi=E_\pi-m_\pi$ and at $|t|$. For an isoscalar nucleus $\pi^+$ and $\pi^-$ scattering are equal, so
  the antineutrino uses the same expression.

### 1.1 Pion–nucleus elastic cross section

The elastic cross section is written as $d\sigma/dt = A_1(T_\pi)\,e^{-b_1(T_\pi)|t|}$ and is built in three pieces.

1. **$T_\pi\le 0.87$ GeV: ¹²C data.** BS09 fitted $A_1$ and $b_1$ to the elastic cross sections that follow
   from the phase-shift analyses of the measured $\pi^\pm$–¹²C differential and total cross sections
   for $T_\pi$ = 30–870 MeV (O. Dumbrajs *et al.*, PRC 29, 581 (1984); H.G. Schlaile, PRC 55, 2584 (1997)).
   The table is from BS09 Fig. 2 (right):

   | $T_\pi$ [GeV] | 0.076 | 0.080 | 0.100 | 0.148 | 0.162 | 0.226 | 0.486 | 0.584 | 0.662 | 0.776 | 0.870 |
   |---|---|---|---|---|---|---|---|---|---|---|---|
   | $A_1$ [mb/GeV²] | 11600 | 14700 | 18300 | 21300 | 22400 | 16400 | 5730 | 4610 | 4570 | 4930 | 5140 |
   | $b_1$ [GeV⁻²] | 116.0 | 109.0 | 89.8 | 91.0 | 89.2 | 80.8 | 54.6 | 55.2 | 58.4 | 60.5 | 62.2 |

   - Between nodes, $A_1$ and $b_1$ are interpolated linearly, as in BS09.
   - Below 76 MeV, where the exponential form fails (BS09), $A_1$ falls linearly to zero at $T_\pi=0$ with
     $b_1$ fixed. This is a negligible part of the cross section.
   - In the Δ region, $\sigma_\text{el}$ from the data peaks at about 250 mb. The Rein–Sehgal absorptive
     model predicts about 1100 mb there, which is the reason for using the data at low energy.
2. **$T_\pi\ge 1.7$ GeV: Rein–Sehgal hadronic model** (D. Rein and L.M. Sehgal, Nucl. Phys. B 223, 29 (1983); BS09 Eqs. 10–15):
   $$\frac{d\sigma}{dt}=A^2\,\frac{\bar\sigma_\text{tot}^2}{16\pi}\,e^{-b|t|}\,F_\text{abs},\quad b=\tfrac13R_0^2A^{2/3},\quad
   F_\text{abs}=\exp\!\Big(-\frac{9A^{1/3}}{16\pi R_0^2}\,\bar\sigma_\text{inel}\Big)$$
   - $R_0 = 1.057$ fm (BS09), so $b = 50.1$ GeV⁻² for ¹²C.
   - $\bar\sigma_\text{tot}=\tfrac12[\sigma(\pi^+p)+\sigma(\pi^-p)]$ comes from the PDG Regge + $\ln^2 s$ fit
     ($Z=20.86$ mb, $Y_1=19.24$ mb, $\eta_1=0.4473$, $B=\pi(\hbar c)^2/M^2$ with $M=2.1206$ GeV; the
     $Y_2$ term cancels in the average). It gives 26.0 mb at $p=10$ GeV and about 24 mb at 100 GeV.
   - $\bar\sigma_\text{inel}=\bar\sigma_\text{tot}-\bar\sigma_\text{el}$. Here
     $\bar\sigma_\text{el}=\bar\sigma_\text{tot}^2/(16\pi B_\text{el})$ is the forward diffraction peak,
     with $B_\text{el}=5.5+0.5\ln(s/\text{GeV}^2)$ GeV⁻² (A. Lebiedowicz and A. Szczurek, PRD 81, 036003
     (2010), arXiv:0912.0190; πp slopes of 7–10 GeV⁻² at 3–200 GeV/c). This gives
     $\bar\sigma_\text{el}$ = 4.9 mb at 10 GeV and 6.7 mb at 3 GeV, and $F_\text{abs}\approx0.45$.
   - The resulting $\sigma_\text{el}(\pi\,^{12}\text{C})$ is 60 mb at $T_\pi = 1.7$ GeV, 50 mb at 5 GeV and 43 mb at 20 GeV.
3. **$0.87<T_\pi<1.7$ GeV: interpolation.** $A_1$ and $b_1$ are interpolated linearly between the last
   data node and the RS value at 1.7 GeV. BS09 extend the last data line to 1.7 GeV and quote ±6% on
   σ(2 GeV) for a 30% change there. BS09 note that the RS and data curves meet near $|\mathbf p_\pi|\approx1$
   GeV at $\sigma_\text{el}\approx80$ mb. Our RS value at $T_\pi$ = 0.87 GeV is 68 mb, because the smooth
   Regge fit has no resonance structure there; the data value is 83 mb.

**Other nuclei.** These are used for validation only; the loop target is CH.
- In the data region the ¹²C table is scaled as $A_1\propto A^{4/3}$, $b_1\propto A^{2/3}$, i.e.
  $\sigma_\text{el}\propto A^{2/3}$ (BS09 section 5).
- Above 1.7 GeV the RS formula is used with the actual $A$.

### 1.2 Parameters

None of these parameters is fitted to the corpus.

| Parameter | Value | Source |
|---|---|---|
| $f_\pi$ | 0.1307 GeV | BS09 |
| $m_A$ | 0.95 GeV | BS09 Eq. 5 |
| $\cos\theta_C$ | $|V_{ud}|$ = 0.97367 | PDG 2024 |
| $R_0$ | 1.057 fm | BS09 Eq. 12 |
| data/RS boundaries | 0.87 / 1.7 GeV | BS09 table range / BS09 section 5 |
| $m_\pi$ | $m_{\pi^\pm}$ | |

YAML keys (`channels:` entry, all optional): `f_pi`, `mA`, `lepton_mass_correction` (true),
`elastic: {R0, table_max_T, rs_min_T}`, `proposal_slope` (0.8).

## 2. Kinematics and sampling

The variables are $(\nu, Q^2, |t|, \phi_\ell, \phi_\pi)$. The weight is
$d\sigma/(dQ^2\,dy\,dt)\,/\,(E\,g_\nu g_{Q^2} g_t)$ in $10^{-38}$ cm² per nucleus (`Multiplicity` = 1).

1. **$\nu$.** Log-uniform in $[m_\pi, E-m_\mu]$, which absorbs the $1/y$ of the integrand.
2. **$Q^2$.** Sampled between the exact lepton limits $2E(E'\mp p')-m_\mu^2$ from
   $g\propto (Q^2+m_A^2)^{-2}$, the $G_A^2$ shape. The low-$Q^2$ region where the pion-pole term acts is
   fully covered. The correction uses BS09's $Q^2_\text{min}=m^2y/(1-y)$; $(Q^2-Q^2_\text{min})$ is clipped at 0.
3. **$|t|$.** $|t| = 0.9\,t_0 + \text{Exp}(b_0)$. Here $t_0=\big[(Q^2+m_\pi^2)/(|\mathbf q|+p_\pi^{(0)})\big]^2$
   is the recoil-free collinear minimum and $b_0 = 0.8\,b_\text{RS}(A)$ is smaller than every $b_1$.
4. **Exact elastic recoil of the ground-state nucleus.** $E_\pi=\nu-|t|/2M_A$. The pion angle to
   $\mathbf q$ follows from $t=-Q^2+m_\pi^2-2(\nu E_\pi-|\mathbf q|p_\pi\cos\theta_{\pi q})$, and trials with
   $|\cos\theta_{\pi q}|>1$ are rejected, so the true $t_\text{min}$ (with recoil) is respected. The
   recoil $P_A' = q+P_A-p_\pi$ then has $P_A'^2 = M_A^2$ exactly.
5. **Azimuths.** The lepton azimuth and the pion azimuth around $\mathbf q$ are uniform. The PCAC
   formula carries no azimuthal dependence.

**Event record.**
- The lepton and the pion are written with status 1. The pion is **not** handed to the cascade:
  distortion and absorption of the outgoing pion are already contained in the measured elastic π–A
  amplitude (the data are elastic, i.e. pion survival).
- The recoiling nucleus is written as the remnant (status 23). The cascade, which has no hadron to
  transport, re-emits it as the status-25 remnant with the balancing four-momentum, which is identical
  to the recoil. It is therefore never a status-1 particle, and the measurement bridge does not count it
  in $E_\text{avail}$ or in the proton/photon vetoes.
- No struck nucleon is recorded (`hit_pdg` = 0). The NuHepMC writer attaches the recoiling nucleus to
  the primary vertex.
- `W` is filled in the free-nucleon convention $W^2=M^2+2M\nu-Q^2$, like the other channels.
  `genpred.py` maps process 100 to the GENIE coherent type (4).
- Trial efficiency is about 2% (8 s for 2M coherent-only events). In the full loop configuration the
  total generation time is unchanged (26 s).

## 3. Validation

**Implementation.**
- `tests/test_coherent.cc` checks the table, the interpolation continuity, the RS model, the PDG πN
  fit, the $A^{2/3}$ scaling, the lepton-mass factor ($m_\ell\to0$ limit and Adler suppression at
  $Q^2_\text{min}$), and events. For events it requires four-momentum conservation to $10^{-9}$, a
  recoil mass of $M_A$, exactly μ + π± with status 1, no nucleons, process 100, $\langle|t|\rangle<0.05$
  GeV², no channel on hydrogen, and σ rising with E.
- In generated events, the slope of $|t|-t_\text{min}$ reproduces the input $b_1$: 91.0 / 89.7 / 80.4 /
  54.3 / 62.4 GeV⁻² against 89.8 / 90.5 / 80.8 / 54.6 / 62.2 at $T_\pi$ = 0.10 / 0.155 / 0.225 / 0.485 /
  0.865 GeV, and 50.0 against 50.1 above 3 GeV.
- σ(ν ¹²C) = 5.3×10⁻⁴⁰ cm² at 1 GeV and 13.2×10⁻⁴⁰ at 2 GeV, which reproduces the BS09 calculation
  (their Fig. 3b, carbon-data curve). This is the same formula, so it checks the implementation, not the physics.

**σ(E) per ¹²C nucleus, ν_μ CC [10⁻⁴⁰ cm²]:** 1.5 (0.6 GeV), 3.9 (0.85), 5.3 (1), 13.2 (2), 18.9 (3), 24.0 (4), 32.1
(6), 44.0 (10), 62.8 (20), 69.2 (25). The ν̄ values are equal within 1% above 1 GeV. The A-dependence at
6 GeV is O/C = 1.15, Ar/C = 1.82, Fe/C = 2.12, i.e. between $A^{1/3}$ and $A^{2/3}$, because heavier
nuclei absorb more through $F_\text{abs}$.

**External data, non-MINERvA.**
- NC is obtained from the same code as $\sigma(\nu_e\,\text{CC}, m_e\approx0)/(2\cos^2\theta_C)$, which is BS09 Eq. 7.
- Other nuclei use the model with the actual $A$, interpolated in $\log A$ between C, O, Ar and Fe.
- Each experiment is evaluated at its mean energy (no flux folding), so the comparison is approximate.

| Experiment | Process, target, ⟨E⟩ | Data | Model |
|---|---|---|---|
| T2K (arXiv:1604.04406) | CC ν, C, 0.85 GeV (E<1.5) | 3.9 ± 1.0 (+1.5 −1.4) ×10⁻⁴⁰ | 3.9×10⁻⁴⁰ |
| ArgoNeuT (arXiv:1408.0598) | CC ν̄, Ar, 3.6 GeV | 5.5 (+2.6 −2.1) ×10⁻³⁹ | 3.6×10⁻³⁹ |
| ArgoNeuT | CC ν, Ar, 9.6 GeV | 2.6 (+1.2 −1.0) ×10⁻³⁸ | 0.87×10⁻³⁸ |
| NOvA (via arXiv:2205.06642 Table 2) | NC, A=13.8, 2.7 GeV | 1.4 ± 0.2 ×10⁻³⁹ | 1.1×10⁻³⁹ |
| MINOS (idem) | NC, 80% Fe + 20% C, 4.9 GeV | 3.26 ± 0.21 ×10⁻³⁹ | 3.0×10⁻³⁹ |
| SKAT (Z. Phys. C 31, 203 (1986)) | CC ν, freon A=30, 7 GeV | 10.6 ± 1.6 ×10⁻³⁹ | 5.8×10⁻³⁹ |
| SKAT | NC, A=30, 7 GeV | 5.2 ± 1.9 ×10⁻³⁹ | 3.3×10⁻³⁹ |
| Aachen–Padova (PLB 125, 230 (1983)) | NC, Al, 2 GeV | 29 ± 10 ×10⁻⁴⁰ | 10.5×10⁻⁴⁰ |
| NOMAD (arXiv:0910.0062) | NC, A=12.8, 24.8 GeV ($E_\pi>0.5$) | 72.6 ± 10.7 ×10⁻⁴⁰ | 39×10⁻⁴⁰ |
| CHARM (PLB 157, 469 (1985)) | NC, A=20, 31 GeV | 96 ± 42 ×10⁻⁴⁰ | 59×10⁻⁴⁰ |
| CHARM II (PLB 313, 267 (1993)) | CC ν, glass A≈20.7, 23.8 GeV | 51 ± 10 ×10⁻⁴⁰ | 100×10⁻⁴⁰ |

- **Below ~5 GeV** (T2K, MINOS, NOvA, ArgoNeuT ν̄, the carbon-data region of the model) the model agrees
  within 1–1.5σ.
- **At 7–25 GeV** the model is below SKAT (×0.55), NOMAD (×0.54) and ArgoNeuT ν (×0.33, 1.7σ). It is
  above CHARM II (×2). These data are mutually inconsistent at the factor-2 level, and the model lies
  inside their spread, on the low side of the majority.
- The high-energy normalisation is set by the RS absorptive optical model.
  - With $F_\text{abs}\approx0.45$, it gives $\sigma_\text{el}(\pi C)\approx45$–50 mb.
  - Measured high-energy π–C total and absorption cross sections leave room for a larger elastic part.
    This has not been checked quantitatively here.
  - This is the main physics uncertainty of the channel at MINERvA energies.
- MINERvA coherent measurements (arXiv:1409.3835, 1711.01178) are deliberately not used here.

**MINERvA ME flux average (CH target).**
- $\langle\sigma_\text{coh}\rangle$ = 0.303×10⁻³⁸ cm² per CH = **2.33×10⁻⁴⁰ cm² per nucleon**, 0.61% of the
  CC total (49.65×10⁻³⁸ per CH).
- In coherent events ⟨E_ν⟩ = 7.0 GeV, ⟨Q²⟩ = 0.18 GeV² (51% below 0.1 GeV²), ⟨|t|⟩ = 0.024 GeV² (99.3%
  below 0.1), ⟨ν/E⟩ = 0.25 and ⟨E_π⟩ = 1.6 GeV, with E_π/ν = 0.999.
- The median θ_π is 14° (90% below 35°), the median θ_μ is 3.2° (90% below 9°), and the median $p_T^\mu$ is 0.26 GeV.
- 23% of coherent pions have $T_\pi<350$ MeV, mostly at 125–350 MeV.

## 4. Limitations

- **Transverse response.** $\sigma_{L,R}$ is neglected (PCAC at $Q^2\to0$ extrapolated with $G_A$),
  together with any azimuthal correlation between the lepton and pion planes.
- **High-energy hadronic model.** Above 1.7 GeV the π–A elastic cross section comes from the RS absorptive
  optical model (one radius parameter, a classical attenuation factor and a smooth πN fit). The
  0.87–1.7 GeV region is an interpolation, and the low-energy table is for ¹²C only (other nuclei are
  scaled by $A^{2/3}$). Glauber-model π–A amplitudes (e.g. Bel'kov and Kopeliovich 1987; arXiv:2205.06642)
  would be the next refinement.
- **No NC coherent channel** (process 150). The loop is CC only.
- **No diffractive scattering on hydrogen.** At high W this is a separate (small) process.
- **Isoscalar target assumed.** π⁺ and π⁻ elastic cross sections are taken as equal (exact for ¹²C).
