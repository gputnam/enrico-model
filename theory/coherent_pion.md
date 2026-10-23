# Charged-current coherent pion production (`CCCOH`, process 100)

Stage 1 `CCCOHChannel`, stage 2 `GroundStateNucleus`, stage 3 `CoherentPionKinematics` (+ `PionNucleusElastic`, `GlauberPionNucleus`),
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

Since model v8 the elastic cross section is built in two pieces: the measured ¹²C data up to
$T_\pi$ = 0.87 GeV and a Glauber multiple-scattering amplitude above. The v7 construction
(Rein–Sehgal optical model above 1.7 GeV, linear bridge in between) is kept as the option
`elastic: {model: rein_sehgal}` and is described at the end of this section.

1. **$T_\pi\le 0.87$ GeV: ¹²C data.** $d\sigma/dt = A_1(T_\pi)\,e^{-b_1(T_\pi)|t|}$. BS09 fitted $A_1$ and $b_1$ to the elastic cross sections that follow
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
2. **$T_\pi > 0.87$ GeV: Glauber amplitude** (`kinematics/GlauberPionNucleus`). R.J. Glauber, in *Lectures in
   Theoretical Physics* Vol. 1 (1959) 315; W. Czyż and L.C. Maximon, Ann. Phys. 52, 59 (1969); the form used for
   hadron–nucleus total cross sections in U. Dersch *et al.* (SELEX), Nucl. Phys. B 579, 277 (2000)
   [hep-ex/9910052], Eq. 24; applied to coherent pion production by H. Sogarwal and P. Shukla, Nucl. Phys. A
   1024, 122494 (2022) [arXiv:2205.06642]. With independent nucleons,
   $$\Gamma(b) = 1-\Big[1-\tfrac12\,\sigma_{\pi N}(1-i\rho)\,T_\beta(b)\Big]^A,\qquad
   T_\beta(b)=\frac{1}{2\pi}\int_0^\infty q\,dq\,J_0(qb)\,S(q)\,e^{-\beta q^2/2},$$
   $$\frac{d\sigma}{dt}=\pi\,\Big|\int_0^\infty b\,db\,J_0(qb)\,\Gamma(b)\Big|^2 e^{q^2 b_{HO}^2/2A},\quad q=\sqrt{|t|},$$
   $$\sigma_\text{tot}=2\!\int\! d^2b\,\text{Re}\,\Gamma,\quad \sigma_\text{el}=\!\int\! d^2b\,|\Gamma|^2,\quad
   \sigma_\text{prod}=\!\int\! d^2b\,\{1-[1-\sigma^\text{in}_{\pi N}T_\beta(b)]^A\}.$$
   The inputs are all measured and none is fitted to neutrino data:
   - **Ground state.** $S(q)$ is the single-nucleon (point) form factor. For ¹²C it is the harmonic-oscillator charge
     density from elastic electron scattering, $\rho_\text{ch}\propto(1+\alpha r^2/a^2)e^{-r^2/a^2}$ with $a$ = 1.687 fm,
     $\alpha$ = 1.067 (H. de Vries, C.W. de Jager, C. de Vries, At. Data Nucl. Data Tables 36, 495 (1987)), so
     $F_\text{ch}(q)=[1-\alpha q^2a^2/(2(2+3\alpha))]e^{-q^2a^2/4}$. The nucleon charge radii are unfolded,
     $S = F_\text{ch}\,e^{q^2(r_p^2+r_n^2)/6}$ with $r_p^2+r_n^2$ = 0.7071 − 0.1161 fm². The fixed-centre shell-model density
     needed by the product wave function carries the extra factor $e^{-q^2b_{HO}^2/4A}$, removed again from the amplitude
     (Czyż–Maximon c.m. correction), with $b_{HO}^2=a^2-\frac23(r_p^2+r_n^2)$, i.e. $b_{HO}$ = 1.57 fm. Then
     $S(q)e^{-\beta q^2/2}=(1-cq^2)e^{-\gamma q^2}$ and $T_\beta(b)$ is analytic.
     Other nuclei (validation only): ¹⁶O HO ($a$ = 1.833 fm, $\alpha$ = 1.544); ²⁷Al, ²⁸Si, ⁴⁰Ar, ⁵⁶Fe two-parameter
     Fermi from de Vries 1987; other $A$ a 2pF with the global rms-radius systematics and $d$ = 0.537 fm. The 2pF
     densities are used as point densities.
   - **πN amplitude.** For the isoscalar ¹²C, $f_{\pi N}=\frac{k\sigma_{\pi N}}{4\pi}(i+\rho)e^{-\beta q^2/2}$ with
     $\sigma_{\pi N}=\frac12[\sigma(\pi^+p)+\sigma(\pi^-p)]$ from the PDG data compilation (RPP 2022 files `pipp_total`,
     `pimp_total`; error-weighted means in 0.05 GeV/c bins below 2.5 GeV/c, coarser above). It includes the N*
     structure: 36–37 mb at $p_\pi$ = 1–1.5 GeV/c, where the smooth Regge fit gives 29–31 mb. The isoscalar
     cross sections are Fermi-averaged with a Gaussian in $s$ of rms $2p_\pi p_F/\sqrt5$, $p_F$ = 0.221 GeV.
   - **Re/Im.** $\rho$ follows from derivative analyticity applied to the PDG Regge + $\ln^2s$ fit:
     $\text{Re}=-\cot(\pi\alpha_R/2)\,Y_1s^{-\eta_1}+\pi B\ln(s/s_M)$ with $\alpha_R=1-\eta_1$. It gives
     $\rho$ = −0.34 / −0.30 / −0.25 / −0.14 / −0.05 at $p$ = 1 / 1.8 / 3.1 / 10 / 31 GeV/c. The measured values are
     π⁻p −0.07…−0.17 at 1.9–6 GeV/c, π⁺p −0.39 at 3.6 GeV/c, and both −0.1…0 at 30–50 GeV/c (PDG `.reim` files).
   - **πN slope.** $\beta=\sigma_\text{tot}^2(1+\rho^2)/(16\pi\sigma_\text{el})$ is taken from the measured isoscalar elastic
     cross section (PDG `pipp_elastic`, `pimp_elastic`, same binning and smearing). This gives 4.8 GeV⁻² at 1 GeV/c,
     7.3 at 10 GeV/c and 8.6 at 60 GeV/c. $\sigma^\text{in}_{\pi N}=\sigma_\text{tot}-\sigma_\text{el}$.
   - **Numerics.** $b$ grid 0–16 fm (800 points). For each $A$, on first use, the complex amplitude is tabulated on
     200 $\ln T$ nodes (0.6–250 GeV) × 321 nodes in $q=\sqrt{|t|}$ (0–1 GeV). This takes 0.13 s for ¹²C and is
     thread safe. Bilinear interpolation of the amplitude, which is smooth through the diffraction minimum,
     reproduces the direct evaluation to < 0.1%.
3. **Matching.** Over $0.87<T_\pi<0.97$ GeV, $d\sigma/dt$ is the linear mixture
   $(1-w)\,A_1e^{-b_1|t|}|_{0.87}+w\,(d\sigma/dt)_\text{Glauber}$, $w=(T_\pi-0.87)/0.1$. Without any adjustment the
   Glauber amplitude at 0.87 GeV already gives $\sigma_\text{el}$ = 83.3 mb and a forward slope of 62.5 GeV⁻²,
   against 82.6 mb and $b_1$ = 62.2 GeV⁻² of the data fit. The blend only removes a sub-percent step.
   Below 0.8 GeV the Glauber model falls below the data (60 against 82 mb at 0.6 GeV): the eikonal/impulse
   approximation fails near the Δ, which is why the data are kept there.

**Hadronic validation of the Glauber piece (¹²C, no free parameter).**

| | Glauber | Data |
|---|---|---|
| σ_tot(π⁻C), 591 GeV/c | 233 mb | 234.1 ± 3.5 mb (SELEX, hep-ex/9910052) |
| σ_tot(π⁺C), 175 / 70 GeV/c | 233 / 231 mb | 244.9 ± 14 / 244.7 ± 40 mb (A. Schiz *et al.*, PRD 21, 3010 (1980)) |
| σ_prod(π⁺C), 60 / 31 GeV/c | 173.5 / 175.2 mb | 171.6 ± 4.6 / 158.3 ± 3.3 mb (NA61/SHINE, PRD 98, 052001 (2018); HEPData ins1672945) |
| σ_inel = σ_tot − σ_el, 60 / 31 GeV/c | 190.4 / 192.1 mb | 188.2 ± 4.7 / 177.0 ± 3.2 mb (NA61) |
| dσ/dt(π±C), 70 and 175 GeV/c, 0.005 < \|t\| < 0.035 GeV² | | data/model 0.80–1.11 (mean ≈ 0.98) (Schiz *et al.*, HEPData ins147258 Tables 4, 22, 49, 67) |
| local slope 0.005–0.04 GeV², 70/175 GeV/c | 66.8 GeV⁻² | 56.6–61.8 GeV⁻² (data include unresolved quasi-elastic scattering, which also gives the data excess above \|t\| ≈ 0.04) |

- The Glauber model includes no inelastic (Gribov) screening, which lowers σ_tot by a few % at high energy
  (SELEX). This is consistent with the 0–8% overshoot of σ_inel/σ_prod at 31 GeV/c.
- At $|t|<0.004$ GeV² the measured dσ/dt is dominated by Coulomb scattering.

**Before → after** (v7 RS bridge → v8 Glauber, ¹²C; σ_el integrated over |t| < 1 GeV²):

| $T_\pi$ [GeV] | 1.0 | 1.7 | 3.0 | 10 | 20 |
|---|---|---|---|---|---|
| σ_el [mb] | 79.7 → 85.8 | 60.4 → 76.2 | 54.6 → 65.2 | 45.9 → 49.1 | 43.2 → 44.3 |
| dσ/dt(0) [mb/GeV²] | 4810 → 5520 | 3030 → 4930 | 2740 → 4250 | 2300 → 3160 | 2170 → 2860 |
| slope 0–0.02 / 0.02–0.06 GeV² [GeV⁻²] | 60.3 → 62.5/78.8 | 50.1 → 62.8/78.6 | 50.1 → 63.3/78.4 | 50.1 → 62.4/75.8 | 50.1 → 62.7/75.8 |

The ¹²C diffraction pattern is not exponential. The local slope steepens from about 63 GeV⁻² at $t\to0$ to about 78 GeV⁻²
at 0.02–0.06 GeV², and the first minimum lies near $|t|\approx$ 0.1 GeV². The RS slope $b=R_0^2A^{2/3}/3$ = 50 GeV⁻²
is too small for carbon: the BS09 data fit itself has $b_1$ = 62 GeV⁻² at 0.87 GeV, and so do the high-energy
π–C data. The integrated σ_el changes little above 5 GeV, because the RS absorptive model already reproduced
σ_el(πC) ≈ 41–46 mb there.

**Option `model: rein_sehgal` (model v7).**
- **$T_\pi\ge 1.7$ GeV: Rein–Sehgal hadronic model** (D. Rein and L.M. Sehgal, Nucl. Phys. B 223, 29 (1983); BS09 Eqs. 10–15):
  $$\frac{d\sigma}{dt}=A^2\,\frac{\bar\sigma_\text{tot}^2}{16\pi}\,e^{-b|t|}\,F_\text{abs},\quad b=\tfrac13R_0^2A^{2/3},\quad
  F_\text{abs}=\exp\!\Big(-\frac{9A^{1/3}}{16\pi R_0^2}\,\bar\sigma_\text{inel}\Big)$$
  - $R_0 = 1.057$ fm (BS09), so $b = 50.1$ GeV⁻² for ¹²C.
  - $\bar\sigma_\text{tot}$ comes from the PDG Regge + $\ln^2 s$ fit ($Z=20.86$ mb, $Y_1=19.24$ mb, $\eta_1=0.4473$,
    $M=2.1206$ GeV). It gives 26.0 mb at $p=10$ GeV.
  - $\bar\sigma_\text{inel}=\bar\sigma_\text{tot}-\bar\sigma_\text{tot}^2/(16\pi B_\text{el})$, with
    $B_\text{el}=5.5+0.5\ln(s/\text{GeV}^2)$ GeV⁻² (A. Lebiedowicz and A. Szczurek, PRD 81, 036003 (2010)). This gives $F_\text{abs}\approx0.45$.
- **$0.87<T_\pi<1.7$ GeV:** $A_1$ and $b_1$ are interpolated linearly between the last data node and the RS value at 1.7 GeV.

**Other nuclei.** These are used for validation only; the loop target is CH.
- In the data region the ¹²C table is scaled as $A_1\propto A^{4/3}$, $b_1\propto A^{2/3}$, i.e.
  $\sigma_\text{el}\propto A^{2/3}$ (BS09 section 5).
- Above 0.87 GeV the Glauber amplitude (or the RS formula) is used with the actual $A$ and density.

### 1.2 Parameters

None of these parameters is fitted to the corpus.

| Parameter | Value | Source |
|---|---|---|
| $f_\pi$ | 0.1307 GeV | BS09 |
| $m_A$ | 0.95 GeV | BS09 Eq. 5 |
| $\cos\theta_C$ | $|V_{ud}|$ = 0.97367 | PDG 2024 |
| $R_0$ | 1.057 fm | BS09 Eq. 12 |
| data/Glauber boundary | 0.87 GeV (+0.1 GeV blend) | BS09 table range |
| ¹²C HO charge density | $a$ = 1.687 fm, $\alpha$ = 1.067 | de Vries *et al.* 1987 (electron scattering) |
| $r_p^2$, $r_n^2$ | 0.7071, −0.1161 fm² | CODATA 2018 / PDG |
| $\sigma_{\pi N}$, $\sigma^\text{el}_{\pi N}$ | PDG data tables, Fermi-averaged ($p_F$ = 0.221 GeV) | PDG RPP 2022 |
| $\rho_{\pi N}$ | Regge-fit analyticity | PDG Regge fit |
| RS option: $R_0$, data/RS boundaries | 1.057 fm; 0.87 / 1.7 GeV | BS09 |
| $m_\pi$ | $m_{\pi^\pm}$ | |

YAML keys (`channels:` entry, all optional): `f_pi`, `mA`, `lepton_mass_correction` (true),
`elastic: {model (glauber | rein_sehgal), table_max_T (0.87), blend_width (0.1), R0, rs_min_T}`, `proposal_slope` (0.8).

## 2. Kinematics and sampling

The variables are $(\nu, Q^2, |t|, \phi_\ell, \phi_\pi)$. The weight is
$d\sigma/(dQ^2\,dy\,dt)\,/\,(E\,g_\nu g_{Q^2} g_t)$ in $10^{-38}$ cm² per nucleus (`Multiplicity` = 1).

1. **$\nu$.** Log-uniform in $[m_\pi, E-m_\mu]$, which absorbs the $1/y$ of the integrand.
2. **$Q^2$.** Sampled between the exact lepton limits $2E(E'\mp p')-m_\mu^2$ from
   $g\propto (Q^2+m_A^2)^{-2}$, the $G_A^2$ shape. The low-$Q^2$ region where the pion-pole term acts is
   fully covered. The correction uses BS09's $Q^2_\text{min}=m^2y/(1-y)$; $(Q^2-Q^2_\text{min})$ is clipped at 0.
3. **$|t|$.** $|t| = 0.9\,t_0 + \text{Exp}(b_0)$. Here $t_0=\big[(Q^2+m_\pi^2)/(|\mathbf q|+p_\pi^{(0)})\big]^2$
   is the recoil-free collinear minimum and $b_0 = 0.8\,b_\text{RS}(A)$ = 40 GeV⁻² (¹²C) is smaller than every $b_1$
   and every local Glauber slope.
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
  0.865 GeV, and 50.0 against 50.1 above 3 GeV with the `rein_sehgal` option.
- σ(ν ¹²C) = 5.3×10⁻⁴⁰ cm² at 1 GeV and 13.2×10⁻⁴⁰ at 2 GeV, which reproduces the BS09 calculation
  (their Fig. 3b, carbon-data curve). This is the same formula, so it checks the implementation, not the physics.

- Glauber (`[coherent]` tests): the analytic HO thickness, σ_tot(π⁻C) and σ_prod(π⁺C) at 60–175 GeV/c, the slope
  and σ_el against the data fit at 0.87 GeV, continuity of the blend, table against direct evaluation (< 0.5%),
  and that the `rein_sehgal` option reproduces v7.

**σ(E) per ¹²C nucleus, ν_μ CC [10⁻⁴⁰ cm²], v7 (RS) → v8 (Glauber):** 1.50 → 1.50 (0.6 GeV), 3.93 → 3.95 (0.85),
5.29 → 5.31 (1), 13.1 → 13.2 (2), 19.1 → 19.5 (3), 24.1 → 25.0 (4), 32.0 → 33.6 (6), 43.7 → 46.1 (10),
62.5 → 65.8 (20), 69.3 → 72.8 (25). The ν̄ values are equal within 1% above 1 GeV.
- The increase is only 1–5%, although σ_el(πC) rises by 7–26% at $T_\pi$ = 1–3 GeV. The steeper
  Glauber slope (63–79 against 50 GeV⁻²) removes strength at finite $|t|$, and in particular near the kinematic
  minimum $t_\text{min}$, so the gain is mostly compensated.
- A-dependence at 6 GeV: O/C = 1.17, Ar/C = 2.02, Fe/C = 2.48. In v7 these were 1.15, 1.82, 2.12; the RS
  attenuation $F_\text{abs}$ over-absorbs heavy nuclei.

**External data, non-MINERvA.**
- NC is obtained from the same code as $\sigma(\nu_e\,\text{CC}, m_e\approx0)/(2\cos^2\theta_C)$, which is BS09 Eq. 7.
- Other nuclei use the model with the actual $A$, interpolated in $\log A$ between C, O, Ar and Fe.
- Each experiment is evaluated at its mean energy (no flux folding), so the comparison is approximate.

| Experiment | Process, target, ⟨E⟩ | Data | v7 (RS) | v8 (Glauber) |
|---|---|---|---|---|
| T2K (arXiv:1604.04406) | CC ν, C, 0.85 GeV (E<1.5) | 3.9 ± 1.0 (+1.5 −1.4) ×10⁻⁴⁰ | 3.9 | 3.9 ×10⁻⁴⁰ |
| ArgoNeuT (arXiv:1408.0598) | CC ν̄, Ar, 3.6 GeV | 5.5 (+2.6 −2.1) ×10⁻³⁹ | 3.6 | 4.1 ×10⁻³⁹ |
| ArgoNeuT | CC ν, Ar, 9.6 GeV | 2.6 (+1.2 −1.0) ×10⁻³⁸ | 0.86 | 1.00 ×10⁻³⁸ |
| NOvA (via arXiv:2205.06642 Table 2) | NC, A=13.8, 2.7 GeV | 1.4 ± 0.2 ×10⁻³⁹ | 1.10 | 1.13 ×10⁻³⁹ |
| MINOS (idem) | NC, 80% Fe + 20% C, 4.9 GeV | 3.26 ± 0.21 ×10⁻³⁹ | 3.0 | 3.6 ×10⁻³⁹ |
| SKAT (Z. Phys. C 31, 203 (1986)) | CC ν, freon A=30, 7 GeV | 10.6 ± 1.6 ×10⁻³⁹ | 5.9 | 6.7 ×10⁻³⁹ |
| SKAT | NC, A=30, 7 GeV | 5.2 ± 1.9 ×10⁻³⁹ | 3.4 | 3.9 ×10⁻³⁹ |
| Aachen–Padova (PLB 125, 230 (1983)) | NC, Al, 2 GeV | 29 ± 10 ×10⁻⁴⁰ | 10.5 | 11.2 ×10⁻⁴⁰ |
| NOMAD (arXiv:0910.0062) | NC, A=12.8, 24.8 GeV ($E_\pi>0.5$) | 72.6 ± 10.7 ×10⁻⁴⁰ | 39 | 41 ×10⁻⁴⁰ |
| CHARM (PLB 157, 469 (1985)) | NC, A=20, 31 GeV | 96 ± 42 ×10⁻⁴⁰ | 63 | 69 ×10⁻⁴⁰ |
| CHARM II (PLB 313, 267 (1993)) | CC ν, glass A≈20.7, 23.8 GeV | 51 ± 10 ×10⁻⁴⁰ | 104 | 114 ×10⁻⁴⁰ |

- **Below ~5 GeV** (T2K, NOvA, ArgoNeuT ν̄) the model agrees within 1–1.5σ. MINOS (mostly Fe) moves from
  ×0.92 to ×1.10, because heavy nuclei gain most.
- **At 7–31 GeV** the model rises by 5–16% and stays low: SKAT CC ×0.63 (was 0.56), NOMAD ×0.57 (0.54),
  ArgoNeuT ν ×0.39 (0.33), CHARM ×0.71 (0.65). CHARM II moves further away (×2.2, was ×2.05). The high-energy
  data still disagree with each other at the factor-2 level.
- **The hadronic input is no longer the open question.** σ_el(πC) above 5 GeV is fixed by the
  measured π–C total, production and elastic dσ/dt data (table in section 1.1), and the Glauber amplitude
  reproduces them without free parameters. The remaining factor ≈2 deficit against SKAT, NOMAD and the MINERvA
  data at $E_\pi>2$ GeV is not a π–A elastic effect. It must come from the PCAC side, i.e. the $Q^2>0$
  extrapolation, the transverse/axial-meson (a₁, ρπ) contributions neglected here, or the coherence-length
  (off-diagonal Glauber) effects of Kopeliovich *et al.*, or from the data definitions.
- Sogarwal and Shukla (arXiv:2205.06642) reach higher high-energy cross sections with a similar Glauber model.
  They do this by using an effective $\alpha_{\pi N}=\text{Re}/\text{Im}$ = 0.9–1.4, fitted to π–nucleus elastic data at
  400–675 MeV, at all energies. The measured forward πN ratio is −0.4…0 at 2–200 GeV/c, so that choice is not
  adopted here.
- MINERvA coherent measurements (arXiv:1409.3835, 1711.01178) are deliberately not used here.

**MINERvA ME flux average (CH target).**
- $\langle\sigma_\text{coh}\rangle$ = 0.317×10⁻³⁸ cm² per CH = **2.44×10⁻⁴⁰ cm² per nucleon** (v7: 0.303, 2.33).
- In coherent events ⟨E_ν⟩ = 7.0 GeV, ⟨Q²⟩ = 0.17 GeV² (v7 0.18), ⟨|t|⟩ = 0.021 GeV² (v7 0.024; the fraction
  above 0.07 GeV² drops from 3.4% to 1.9%), ⟨E_π⟩ = 1.66 GeV (1.62), and the median θ_π is 13.5° (14.3°).

## 4. Limitations

- **Transverse response.** $\sigma_{L,R}$ is neglected (PCAC at $Q^2\to0$ extrapolated with $G_A$),
  together with any azimuthal correlation between the lepton and pion planes.
- **Glauber amplitude.** The following are not included:
  - inelastic (Gribov) screening, a few % on σ_tot at high energy;
  - NN correlations beyond the c.m. correction, and spin-flip πN amplitudes;
  - the longitudinal momentum transfer of the neutrino process: dσ/dt is evaluated at $|t|$ as in BS09,
    without the $e^{iq_Lz}$ phase inside the absorbing nucleus;
  - medium modification of the N* resonances.
  Below 0.8 GeV the eikonal approximation fails, and the data table is used there.
- **ρ_πN in the N* region** comes from the Regge-fit analyticity relation (−0.3 to −0.35 at 1–2 GeV/c). With ρ = 0,
  σ_el(πC) would be 10% lower at 1 GeV and 4% lower at 5 GeV.
- **Low-energy table.** It is for ¹²C only (other nuclei are scaled by $A^{2/3}$). It has no node between 0.226 and
  0.486 GeV, so the falling Δ tail is interpolated linearly.
- **|t| shape.** The Glauber/data slope for ¹²C (62–79 GeV⁻²) is steeper than the RS value (50). The model therefore
  puts less strength at $|t|>0.04$ GeV² than v7. Hadron data support this.
- **No NC coherent channel** (process 150). The loop is CC only.
- **No diffractive scattering on hydrogen.** At high W this is a separate (small) process.
- **Isoscalar target assumed.** π⁺ and π⁻ elastic cross sections are taken as equal (exact for ¹²C).
