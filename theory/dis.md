# Charged-current shallow- and deep-inelastic scattering (CCDIS)

This note covers the `nuflow` model for inclusive inelastic charged-current scattering above the
resonance region, $\nu_\ell N\to\ell^-X$ and $\bar\nu_\ell N\to\ell^+X$ with $W>W_\text{min}$
(default 1.6 GeV, equal to the CCRES `w_max`). It has four parts:

| Stage | Class | Role |
|---|---|---|
| 1 | `channel/CCDISChannel` | Either nucleon can be struck (proton with probability $Z/A$); multiplicity $A$; NuHepMC process id 600 |
| 2 | `initial_state/RelativisticFermiGas`, `FreeNucleon` | Struck-nucleon momentum and removal energy (same as CCQE/CCRES) |
| 3 | `kinematics/DISKinematics`, `kinematics/BodekYang`, `kinematics/GRV98LO` | Lepton kinematics and the weight from nucleon structure functions |
| 3 (after acceptance) | `kinematics/DISHadronization` | Fragmentation of the hadronic system into a nucleon and pions (status 22, then status 1 through stage 4) |

The structure functions are the effective leading-order model of Bodek and Yang [1] (the 2013 version,
v2 of arXiv:1011.6592): GRV98 LO parton densities [2] evaluated at the scaling variable $\xi_w$,
multiplied by low-$Q^2$ $K$ factors. The parameters were fitted by Bodek and Yang to charged-lepton
data only (SLAC, BCDMS, NMC, H1, photoproduction), plus PCAC and high-energy neutrino constraints for
the axial part. Nothing in this channel is fitted to the MINERvA corpus. The code was written from
the papers; the PDF grid is the GRV98 authors' table.

---

## 1. Cross section

### 1.1 Hadronic tensor and structure functions

For a nucleon of four-momentum $P$ (mass $M$) and momentum transfer $q$,
$$
W^{\mu\nu} = -g^{\mu\nu}W_1+\frac{P^\mu P^\nu}{M^2}W_2-\frac{i\epsilon^{\mu\nu\alpha\beta}P_\alpha q_\beta}{2M^2}W_3
+\frac{q^\mu q^\nu}{M^2}W_4+\frac{P^\mu q^\nu+q^\mu P^\nu}{2M^2}W_5 ,
$$
with $\nu=P\cdot q/M$, $x=Q^2/2P\cdot q$ and the dimensionless structure functions
$$
W_1=\frac{F_1}{M},\quad W_2=\frac{F_2}{\nu},\quad W_3=\frac{F_3}{\nu},\quad W_4=\frac{F_4}{x\nu},\quad W_5=\frac{2F_5}{\nu}.
$$
This normalisation of $F_4$, $F_5$ is the one in which the rest-frame cross section has the form of
Kretzer and Reno [6] (and Paschos–Yu), and in which the leading-order Albright–Jarlskog relations read
$F_4=0$, $2xF_5=F_2$. I checked the latter directly: the free-quark tensor
$\propto 2x^2P^\mu P^\nu+x(P^\mu q^\nu+q^\mu P^\nu)-g^{\mu\nu}xP\cdot q$ gives $W_4=0$ and $W_5=W_2/x$.

The lepton tensor is $L^{\mu\nu}=8[k^\mu k'^\nu+k'^\mu k^\nu-g^{\mu\nu}k\cdot k'\mp i\epsilon^{\mu\nu\alpha\beta}k_\alpha k'_\beta]$
(the lepton mass drops out of the trace). The contraction, `DISKinematics::Contract`, is
$$
L_{\mu\nu}W^{\mu\nu}=16(k\cdot k')W_1+\frac{8}{M^2}\big[2(k\cdot P)(k'\cdot P)-M^2k\cdot k'\big]W_2
\pm\frac{8}{M^2}\big[(k\cdot P)(k'\cdot q)-(k\cdot q)(k'\cdot P)\big]W_3
$$
$$
+\frac{8}{M^2}\big[2(k\cdot q)(k'\cdot q)-q^2\,k\cdot k'\big]W_4
+\frac{8}{M^2}\big[(k\cdot P)(k'\cdot q)+(k\cdot q)(k'\cdot P)-(P\cdot q)(k\cdot k')\big]W_5 ,
$$
with the upper sign for neutrinos. For a nucleon at rest,
$$
\frac{d^2\sigma}{dx\,dy}=\frac{G_F^2M\nu}{16\pi E}L_{\mu\nu}W^{\mu\nu}
=\frac{G_F^2ME}{\pi}\Big\{\Big(xy^2+\frac{m^2y}{2EM}\Big)F_1+\Big[1-\frac{m^2}{4E^2}-\Big(1+\frac{Mx}{2E}\Big)y\Big]F_2
\pm\Big[xy\Big(1-\frac y2\Big)-\frac{m^2y}{4EM}\Big]F_3+\frac{m^2(m^2+Q^2)}{4E^2M^2x}F_4-\frac{m^2}{EM}F_5\Big\},
$$
which the test `Tensor contraction reproduces the rest-frame d2sigma/dxdy` verifies to $10^{-9}$
for each $F_i$ separately, with the muon mass, for $\nu$ and $\bar\nu$. The $W$-boson propagator
$(1+Q^2/M_W^2)^{-2}$ is included. CKM factors are in the parton sums (§1.2), not in the prefactor.

### 1.2 Bodek–Yang effective LO structure functions [1]

**Scaling variable** (Eq. 12 of [1]):
$$
\xi_w=\frac{2x\,(Q^2+M_f^2+B)}{Q^2\big[1+\sqrt{1+4M^2x^2/Q^2}\big]+2Ax},
$$
with $M_f=0$ for light final quarks and $M_f=M_c=1.32$ GeV for charm production. $A$ absorbs target-mass,
higher-twist and higher-order effects, $B$ the initial transverse momentum and the final-state effective mass.

**PDFs.** GRV98 LO [2], interpolated from the authors' grid (`data/pdf/grv98lo.dat`: 68 $\ln x$ × 27
$\ln Q^2$ nodes, $0.8<Q^2<10^6$ GeV², linear interpolation in $(\ln x,\ln Q^2)$ after dividing out the
smooth $x$-behaviour, as in the authors' GRV98PA routine). They are evaluated at $(\xi_w,\max(Q^2,0.8))$:
the evolution is frozen below 0.8 GeV² (item 4 of [1], §3). All PDFs are multiplied by $N=1.026$.
The NMC $d/u$ correction (Eqs. 35–37 of [1]) is applied to the valence quarks:
$(d/u)'=d/u+\delta$, $\delta=-0.00817+0.0506\xi+0.0798\xi^2$, with $u_v+d_v$ unchanged. The neutron
uses the isospin-rotated proton PDFs.

**Vector $K$ factors** (second iteration, Eq. 21 and Table 2 of [1]; $G_D=(1+Q^2/0.71)^{-2}$):
$$
K^V_{\text{val},u}=K_{LW}\,[1-G_D^2]\frac{Q^2+C_{v2u}}{Q^2+C_{v1u}},\quad
K^V_{\text{val},d}=K_{LW}\,[1-G_D^2]\frac{Q^2+C_{v2d}}{Q^2+C_{v1d}},\quad
K^V_{\text{sea},f}=\frac{Q^2}{Q^2+C_{\text{sea},f}},\quad K_{LW}=\frac{\nu^2+C_{\text{low-}\nu}}{\nu^2}.
$$
For the neutron, the $u$ and $d$ factors follow the distributions (they swap with them).

**Axial $K$ factors**, "type II" (§11.2 of [1], the default favoured by CCFR/CHORUS):
$K^A_\text{val}=(Q^2+0.018)/(Q^2+0.18)$ and $K^A_\text{sea}=(Q^2+0.018)/(Q^2+0.3)$. They stay finite at
$Q^2=0$ (the PCAC contribution to $F_2$).

**Neutrino structure functions** (Eqs. 22–24 and §8 of [1]). With $\xi q\equiv\xi_wq(\xi_w)$ and
CKM weights $|V_{ij}|^2$:
- $\nu$ (W⁺): quarks $d_v,\bar d$ (as $d$ sea) $\to u$ with $|V_{ud}|^2$; $s\to u$ with $|V_{us}|^2$;
  antiquark $\bar u\to\bar d,\bar s$ with $|V_{ud}|^2+|V_{us}|^2$. Charm: $d\to c$ ($|V_{cd}|^2$),
  $s\to c$ ($|V_{cs}|^2$) at $\xi_w(M_c)$.
- $\bar\nu$ (W⁻): $u_v,\bar u$ (as $u$ sea) with $|V_{ud}|^2+|V_{us}|^2$; antiquarks $\bar d$ ($|V_{ud}|^2$),
  $\bar s$ ($|V_{us}|^2$); charm $\bar d,\bar s\to\bar c$ at $\xi_w(M_c)$.

Then
$$
F_2=\sum K^V\xi q+\sum K^A\xi q\;(\text{quarks and antiquarks}),\qquad
2xF_1=2\,F_2^V\frac{1+4M^2x^2/Q^2}{1+R(x,Q^2)},\qquad
xF_3=2H(x)\Big[\sum_q K^V\xi q-\sum_{\bar q}K^V\xi\bar q\Big],
$$
where $F_2^V$ is the vector sum, $2xF_1^A=2xF_1^V$ (§11 of [1]: the PCAC part is purely longitudinal),
and $H(x)=0.914+0.296x-0.374x^2+0.165x^3$ (Eq. 27). The charm pieces of $2xF_1$ and $xF_3$ carry
$K_\text{charm}=Q^2/(Q^2+M_c^2)$ and $R(\xi_w(M_c),Q^2)$ (Eqs. 30–31). $F_4=0$, $2xF_5=F_2$.

**$R=\sigma_L/\sigma_T$**: R1998 of E143 [3], $R=(R_a+R_b+R_c)/3$ with the Table II coefficients,
for $Q^2\ge0.3$ GeV²; below, $R=3.633\,Q^2/(Q^4+1)\,R_{1998}(x,0.3)$, which vanishes at $Q^2=0$ (§7 of [1]).

**Parameters** (all from [1]; units GeV²):

| $A$ | $B$ | $C_{v1u}$ | $C_{v2u}$ | $C_{v1d}$ | $C_{v2d}$ | $C^u_\text{sea}$ | $C^d_\text{sea}$ | $C^s_\text{sea}$ | $C_{\text{low-}\nu}$ | $N$ | $P^A_\text{val}$, $C^A_\text{val}$ | $P^A_\text{sea}C^A_\text{sea}$, $C^A_\text{sea}$ |
|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.621 | 0.380 | 0.417 | 0.264 | 0.341 | 0.323 | 0.369 | 0.561 | 0.561 | 0.218 | 1.026 | 0.018, 0.18 | 0.018, 0.3 |

CKM: $|V_{ud}|=0.97367$, $|V_{us}|=0.22431$, $|V_{cd}|=0.221$, $|V_{cs}|=0.975$ (PDG 2024).
These are the Table 2 (2013) values, not the earlier 2003/2004 sets ($A=0.419/0.538$, ...).

### 1.3 Bound nucleons

Stage 2 provides $\mathbf p$ and $E=\sqrt{\mathbf p^2+M^2}-\epsilon_b$ (Smith–Moniz). The physical
hadronic mass follows from energy conservation, $W^2=(p_\text{bound}+q)^2$. The structure functions are
taken to be those of a **free nucleon at the same $Q^2$ and $W$**:
$$
x_\text{eff}=\frac{Q^2}{W^2-M^2+Q^2},\qquad \nu_\text{eff}=\frac{W^2-M^2+Q^2}{2M},
$$
and the tensor is built from the on-shell four-momentum $p_\text{on}=(\sqrt{\mathbf p^2+M^2},\mathbf p)$
and the lepton $q$. Binding therefore lowers $W$ at fixed lepton kinematics, and Fermi motion enters
through $p_\text{on}$ and the flux factor. Two alternatives were tried and rejected:
- The QE/RES de Forest choice $\tilde q=q-(\epsilon_b,\mathbf 0)$ inside the tensor violates current
  conservation; the $W_5=W_2/x$ term then picks up a spurious piece $\sim\epsilon_b/(Mx)$. That reduced
  σ on carbon by 15–20% at 40 GeV.
- With $\tilde q$ only in the arguments of $F_i$, $\tilde Q^2=Q^2+2\epsilon_bq_0$ is shifted by
  ~1 GeV² at $q_0\sim20$ GeV. That reduced σ by 7–9%, which is unphysical.

With the adopted prescription, the bound/free ratio is 0.96–0.98 (6 GeV) and 0.985 (40 GeV). With
$\epsilon_b=0$ (Fermi motion only) it is +3%. No EMC/shadowing correction is applied (see §6).

### 1.4 Phase space, flux factor and sampling

As for CCQE/CCRES, the flux factor is $4E_\nu E_N$ (on-shell $E_N$):
$$
d\sigma=\frac{G_F^2M}{32\pi^2E_\nu E_N}\,L_{\mu\nu}W^{\mu\nu}\,\frac{d^3k'}{E'},\qquad
\frac{d^3k'}{E'}=\frac{|\mathbf p^*_\ell|}{2\sqrt s}\,dW^2\,d\Omega^*_\ell ,
$$
in the CM frame of $k+p_\text{bound}$.
- **$W$ proposal**: on $[W_\text{min},\sqrt s-m_\ell]$, half flat in $W$ and half flat in $W^2$.
- **Lepton-angle proposal**: in $x^*=1-\cos\theta^*$, half flat and half $\propto(1+bx^*)^{-2}$ with $b=2k^*p^*_\ell/(1\text{ GeV}^2)$.
- **Weight**: $d\sigma/(dW\,dx^*\,d\phi)$ divided by the proposal density.
- **Efficiency**: on carbon the acceptance $\langle w\rangle/w_\text{max}$ is 0.20 at 6 GeV and 0.13 at 40 GeV.
- **Generation time**: 2M ME-flux events take 13 s; the tables take about 40 s more.

The hadronization (§2) does not change the weight. It runs in `KinematicsGenerator::Finalize`, called
by the generator only for accepted events and before stage 4, so the table building does not pay
for it.

## 2. Hadronization (`DISHadronization`)

The aim is a complete, charge- and baryon-number-conserving final state with the measured mean
charged multiplicity and dispersion. It is not a string model. It was written here and does not
port any generator's scheme.

1. **Charge and baryon.** The hadronic charge is $Q_X=Q_N\pm1$. One baryon is produced, chosen from
   the spectator diquark of valence scattering plus a $u$ or $d$ quark with equal probability (strangeness
   neglected):
   - $\nu p$ ($uu$): always $p$ ($uuu=\Delta^{++}\to p\pi^+$);
   - $\nu n$ and $\bar\nu p$ ($ud$): $p$ or $n$, 50:50;
   - $\bar\nu n$ ($dd$): always $n$.
2. **Charged multiplicity.** The mean is $\langle n_\text{ch}\rangle=a+b\ln W^2$, from deuterium bubble chambers [4,5]:

   | channel | $a$ | $b$ | source |
   |---|---|---|---|
   | $\nu p$ | 0.50 | 1.42 | FNAL E545 [4] |
   | $\nu n$ | −0.20 | 1.42 | FNAL E545 [4] |
   | $\bar\nu p$ | 0.02 | 1.28 | WA25 [5] |
   | $\bar\nu n$ | 0.80 | 0.95 | WA25 [5] |

   The dispersion is $D_\text{ch}=0.35\langle n_\text{ch}\rangle$, the slope measured by E545 [4] (KNO-like). Sampling:
   - write $n_\text{ch}=|Q_X|+2n_\text{pair}$;
   - draw $n_\text{pair}$ by stochastic rounding of a gamma variate with mean $(\langle n_\text{ch}\rangle-|Q_X|)/2$ (floored at 0) and variance $(D_\text{ch}/2)^2-1/6$;
   - raise $n_\text{ch}$ if a charged baryon requires it.
3. **Charged pions.** $n_{\pi^+}-n_{\pi^-}=Q_X-Q_B$, and their sum is $n_\text{ch}$ minus 1 if the baryon is a proton.
4. **Neutral pions.** By isospin, $n_{\pi^0}\sim$ Poisson$\big((n_{\pi^+}+n_{\pi^-})/2\big)$. One $\pi^0$ is forced if no meson was drawn ($W>M+m_\pi$). π⁰s are left undecayed.
5. The species are redrawn if their total mass exceeds $W$.
6. **Kinematics**, in the $X$ rest frame with $z$ along $\mathbf q$:
   - $W<2.3$ GeV, or two bodies: isotropic $n$-body phase space (Raubold–Lynch/GENBOD, unweighted by rejection).
   - $W\ge2.3$ GeV: limited-$p_T$ longitudinal phase space. Gaussian $p_T$ with $\langle p_T^2\rangle=0.25$ GeV² (a typical fragmentation value; a model choice, not fitted), shifted so that $\sum\mathbf p_T=0$. Rapidities are uniform on a plateau, with the baryon at the backward (target-fragmentation) end. They are stretched and shifted so that $\sum E=W$ and $\sum p_z=0$ exactly. The invariant mass $M^2(c)=\sum_{ij}m_{Ti}m_{Tj}\cosh c(y_i-y_j)$ is solved by bisection.

Hadrons are written with status 22 (mother = struck nucleon) and become status 1 in stage 4.

**Achieved multiplicities** (`<n_ch>` incl. protons, 20k decays each):

| | W=1.8 | 2.5 | 3.0 | 4.0 | 6.0 GeV |
|---|---|---|---|---|---|
| $\nu p$: $\langle n_\text{ch}\rangle$ (fit) | 2.16 (2.17) | 3.12 (3.10) | 3.61 (3.62) | 4.44 (4.44) | 5.59 (5.59) |
| $\nu n$ | 1.46 (1.47) | 2.41 (2.40) | 2.92 (2.92) | 3.73 (3.74) | 4.88 (4.89) |
| $D/\langle n\rangle$ ($\nu p$) | 0.25 | 0.36 | 0.36 | 0.35 | 0.35 |
| $\langle n_{\pi^0}\rangle$ ($\nu p$) | 0.57 | 1.05 | 1.31 | 1.71 | 2.29 |

Below $W\approx2$ GeV the charge constraint dominates (for $\nu p$, $n_\text{ch}\ge2$), so the dispersion is smaller.

## 3. Boundary with CCRES/CCNSTAR, single-pion decomposition, process ids

- **Boundary.** CCRES (Δ region, HNV model) and CCNSTAR ($N(1440)$, $N(1520)$, $N(1535)$, [nstar.md](nstar.md)) generate $W<1.6$ GeV; CCDIS generates $W\ge1.6$ GeV. All use the same $W=\sqrt{(p_\text{bound}+q)^2}$, so there is neither a gap nor double counting in $W$.
- **Remaining step.** At the boundary, $d\sigma/dW$ on CH (ME flux, per nucleon) goes from 1.07 (1.5–1.6 GeV) to 1.67 ×10⁻³⁸ cm²/GeV (1.6–1.7 GeV). Without CCNSTAR it went from 0.51–0.53 to 1.66–1.72 at 6 GeV.
- **Process id.** `nuflow` writes NuHepMC process id 600 (CC DIS); the measurement bridge maps it to the DIS interaction type.

### 3.1 Exclusive single-pion part of the continuum

The hadronization of §2 gives a single-pion ($N\pi$) probability $P_{1\pi}(W)$. It is tabulated at construction by Monte Carlo with a fixed seed: 500 nodes in $1.07<W<6.06$ GeV, with a precision of about 0.5%. Its values are:
- $\nu p$: 0.61 at 1.6 GeV, 0.47 at 2 GeV, 0.19 at 3 GeV;
- $\nu n$: 0.77, 0.49 and 0.12.

**Isospin decomposition.** $S(W)=P_{1\pi}\,L_{\mu\nu}W^{\mu\nu}$ is the single-pion strength the continuum implies at a given kinematic point.
- With `single_pion: {isospin: true}`, it is split into isospin components. The pure-$I=3/2$ target ($\nu p$, $\bar\nu n$) has $S_{3/2}=S$.
- On the mixed target ($\nu n$, $\bar\nu p$), $S_{3/2}=\min(S,S_\text{partner}/3)$, where $S_\text{partner}$ is the partner nucleon's single-pion strength at the same kinematics, and $S_{1/2}=S-S_{3/2}$. The $I=1/2$–$3/2$ interference is neglected.

**Scaling and charge states.** Each component is scaled, $\sigma_{1\pi}=R_{3/2}S_{3/2}+R_{1/2}S_{1/2}$, and the Clebsch–Gordan coefficients fix the charge states:
- $I=3/2$: $n\pi^+:p\pi^0=1:2$;
- $I=1/2$: $n\pi^+:p\pi^0=2:1$;
- $\nu p\to p\pi^+$ on the pure target.

**Inclusive strength.** The event is single-pion with probability $\sigma_{1\pi}/(L_{\mu\nu}W^{\mu\nu})$, and otherwise has two or more pions, with the hadronization model's distribution restricted to that class. The weight, and hence the inclusive Bodek–Yang cross section, is unchanged: only the exclusive composition moves.

**Fit of the scales.** $R_{3/2}$ and $R_{1/2}$ are fitted to the ANL/BNL deuterium $p\pi^+$, $n\pi^+$ and $p\pi^0$ cross sections without a W cut, at $E_\nu$ = 1.5–4.5 GeV. The fit model is CCRES + CCNSTAR + scaled CCDIS on free nucleons (details in [nstar.md](nstar.md) §4).
- Result: $R_{3/2}=0.28\pm0.15$ and $R_{1/2}=0.68\pm0.11$, with $\chi^2=11.1/14$ (49.5 at $R=1$).
- The KNO-type hadronization therefore overestimates the single-pion share just above the resonance region, by about 3.5× in $I=3/2$ ($\nu p\to p\pi^+$ at 4.5 GeV would otherwise be 1.28 against 0.995 ± 0.10 measured).
- The loop configuration uses these values. Without `single_pion` the legacy behaviour holds: all multiplicities come from the hadronization, the $W$ dependence is unchanged, and the scales are 1.
- **Not done:** a non-resonant multi-pion continuum below 1.6 GeV, i.e. Bodek–Yang times $1-P_{1\pi}$, which an interrupted earlier attempt tried. Near threshold the KNO $P_{1\pi}$ puts 30–40% of the strength into $N\pi\pi$, so $\nu p\to p\pi^+\pi^0$ comes out 3–10× above the ANL two-pion data (0.014 ± 0.005 at 1.1 GeV, 0.066 ± 0.025 at 1.7 GeV).

## 4. Validation (`tests/test_dis.cc`)

- **GRV98 grid:** at $Q^2=2$, 10 and 100 GeV², the momentum sum $\int x\sum(q+\bar q)+xg$ is 1 within 2%, and $\int u_v=2$, $\int d_v=1$ within 3%.
- **Contraction:** matches the Kretzer–Reno $d^2\sigma/dxdy$ to $10^{-9}$ for each $F_i$, $\nu$ and $\bar\nu$, with the muon mass.
- **Structure-function limits:**
  - at $Q^2=200$, $\xi_w\to x$ and $F_2^{\nu p}=2N(d_v+\bar d+\bar u+s)$;
  - $2xF_1/F_2=1/(1+R)$;
  - as $Q^2\to0$ at fixed $\nu$, $2xF_1/F_2\to0$ while $F_2$ (axial) stays finite;
  - $R_{1998}(0.1,2)\in[0.15,0.4]$.
- **Hadronization:** for ν/ν̄ on p/n and $W$ = 1.7–8 GeV:
  - exact charge, baryon number (one nucleon) and four-momentum conservation;
  - at least one meson;
  - on-shell masses;
  - $\langle n_\text{ch}\rangle$ within 8% of the data fits and $D/\langle n\rangle=0.35\pm0.07$ for $W\ge3$ GeV.
- **Events:** CC events on CH at 6 GeV (ν and ν̄) conserve four-momentum (including the remnant) and charge, and have exactly one lepton, one nucleon, at least one pion and $W\ge1.6$.
- **Free nucleon at 40 GeV:**
  - $\sigma/E=0.424$ ($\nu p$) and 0.844 ($\nu n$), so 0.634 isoscalar for $W>1.6$;
  - independent check by direct $(x,y)$ quadrature: 0.425 and 0.845;
  - carbon per nucleon: 0.624.

**Total $\sigma_\text{CC}/E$ on carbon per nucleon** (QE+RES+DIS, isoscalar), against MINOS on iron [7], which is normalised to the world average 0.675 at 30–50 GeV:

| $E$ [GeV] | 1.08 | 1.92 | 3.41 | 4.13 | 6.06 | 8.89 | 10.8 | 13.1 | 19.2 | 41.3 |
|---|---|---|---|---|---|---|---|---|---|---|
| nuflow (×10⁻³⁸ cm²/GeV) | 0.743 | 0.612 | 0.604 | 0.610 | 0.621 | 0.632 | 0.633 | 0.638 | 0.643 | 0.651 |
| MINOS | | | 0.748±0.061 (3.48) | 0.711±0.033 (4.45) | 0.708±0.032 (5.89) | 0.722±0.045 (7.97) | 0.699±0.043 (10.45) | 0.691±0.028 (13.43) | 0.689±0.016 (19.87) | 0.675±0.018 (38.87) |

The model is 3.5% below the world average at 40 GeV. Bodek–Yang obtain 0.673 for iron, but that
includes their Fe/D nuclear correction, which adds about 3–5% at $0.05<x<0.3$ (anti-shadowing). The
model is about 12% below at 4–8 GeV. Most of that gap is the missing inelastic strength at
$W<1.6$ GeV beyond single-Δ pion production (§6); 2p2h is now the separate `CCMEC` channel ([mec.md](mec.md)).
With `CCMEC` and the second-resonance-region channel `CCNSTAR` ([nstar.md](nstar.md)), the full loop
configuration gives 0.667 (4.47 GeV) and 0.661 (6.50 GeV), 6–7% (≈1.4σ) below MINOS.

**MINERvA ME flux on CH** (2M events), flux-averaged per nucleon, in $10^{-38}$ cm²:

| channel | σ | fraction |
|---|---|---|
| CCQE | 0.407 | 11% |
| CCRES | 0.548 | 15% |
| CCDIS | 2.675 | 74% |
| total | 3.63 | |

DIS events have $\langle Q^2\rangle=2.65$ GeV² and $\langle W\rangle=2.62$ GeV; 25% have $W<2$.

## 5. Configuration

```yaml
- type: CCDIS
  nuclear_model: {type: RFG}
  w_min: 1.6                    # must equal the CCRES w_max
  bodek_yang: {A: 0.621, B: 0.380, N: 1.026, axial: typeII, charm: true, du_correction: true}   # defaults
  hadronization: {dispersion_slope: 0.35, mean_pt2: 0.25, w_lps: 2.3}                           # defaults
  single_pion: {scale_i32: 0.28, scale_i12: 0.68, isospin: true}   # section 3.1; default: none (all 1, no isospin split)
  w_max: 0                      # optional upper W limit (0 = kinematic limit)
  proposal_q2_scale: 1.0
```

## 6. Limitations

- **Resonance–DIS transition.** The cut at $W=1.6$ GeV is still sharp:
  - $d\sigma/dW$ steps up by about 1.6×;
  - the mean π⁰ multiplicity jumps from 0.33 to 0.86;
  - there is no non-resonant multi-pion production below 1.6 GeV (only the $N^*\to N\pi\pi$ decays of CCNSTAR).

  The single-pion scales of §3.1 are constant in $W$ and carry the uncertainties of the fit; the $I=3/2$ scale is weakly constrained (±50%).
- **No nuclear modification of the structure functions** (EMC, anti-shadowing, shadowing) beyond
  Fermi motion and binding in the kinematics. Bodek–Yang's charged-lepton $A/D$ ratios would add about
  +3–5% at $0.05<x<0.3$ and −10% at $x\approx0.6$. Applying them on top of explicit Fermi motion would
  double count.
- **Hadronization:**
  - pions and nucleons only: no strange or charmed hadrons, no η/ω/ρ;
  - charm-production events are fragmented as light-flavour states, so energy and charge are conserved but the flavour is wrong (a few per cent of σ at 6 GeV);
  - no resonance decays in the low-$W$ final state; isotropic phase space below 2.3 GeV;
  - no forward/backward charge correlations beyond the backward baryon;
  - $\langle p_T^2\rangle$ is a single model value.
- **No Pauli blocking** of the outgoing nucleon: negligible at $W>1.6$ GeV.
- **FSI** is applied by the cascade ([fsi_cascade.md](fsi_cascade.md)). The DIS hadrons first travel
  the SKAT formation length $|\vec p|/\mu^2$ ($\mu^2=0.08$ GeV²) before they can interact.
- **Radiative corrections** are not included (Born level, as in [1]).

## References

1. A. Bodek and U.K. Yang, "Axial and vector structure functions for electron- and neutrino-nucleon scattering cross sections at all $Q^2$ using effective leading order parton distribution functions", arXiv:1011.6592 (v2, 2013). Earlier versions: Nucl. Phys. B Proc. Suppl. 112, 70 (2002) [hep-ex/0203009]; hep-ex/0308007; A. Bodek, I. Park, U.K. Yang, Nucl. Phys. B Proc. Suppl. 139, 113 (2005) [hep-ph/0411202].
2. M. Glück, E. Reya, A. Vogt, "Dynamical parton distributions revisited", Eur. Phys. J. C 5, 461 (1998) [hep-ph/9806404]. The grid was transcribed from the DATA-statement version of the authors' GRV98PA routine (R. Engel, 1998).
3. K. Abe et al. (E143), "Measurements of $R=\sigma_L/\sigma_T$ for $0.03<x<0.1$ and fit to world data", Phys. Lett. B 452, 194 (1999) [hep-ex/9808028].
4. D. Zieminska et al. (FNAL E545), "Charged particle multiplicity distributions in νn and νp charged current interactions", Phys. Rev. D 27, 47 (1983).
5. S. Barlag et al. (WA25), Z. Phys. C 11, 283 (1982) [Erratum 14, 281]. Compilation and re-analysis: K.S. Kuzmin and V.A. Naumov, Phys. Rev. C 88, 065501 (2013) [arXiv:1311.4047], Table I.
6. S. Kretzer and M.H. Reno, "Tau neutrino deep inelastic charged current interactions", Phys. Rev. D 66, 113007 (2002); C.H. Albright and C. Jarlskog, Nucl. Phys. B 84, 467 (1975).
7. P. Adamson et al. (MINOS), "Neutrino and antineutrino inclusive charged-current cross section measurements with the MINOS near detector", Phys. Rev. D 81, 072002 (2010) [arXiv:0910.2201], Table III.
8. F. James, "Monte Carlo phase space", CERN 68-15 (1968).
