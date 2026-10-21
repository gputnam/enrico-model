# Charged-current 2p2h excitations by meson-exchange currents (CCMEC)

This note covers the `nuflow` model for charged-current two-particle–two-hole (2p2h) excitations
induced by two-body meson-exchange currents (MEC):

$$
\nu_\ell + (np) \to \ell^- + pp,\qquad \nu_\ell + (nn)\to \ell^- + np,\qquad
\bar\nu_\ell + (pp)\to \ell^+ + np,\qquad \bar\nu_\ell + (np)\to\ell^+ + nn .
$$

| Stage | Class | Role |
|---|---|---|
| 1 | `channel/CCMECChannel` | Pair type drawn per trial (½ each); multiplicity $2N$ (ν) or $2Z$ (ν̄); NuHepMC process id 300 (CC MEC/2p2h range 300–349) |
| 2 | `initial_state/RelativisticFermiGas` | Two hole momenta, each uniform in its Fermi sphere, each with the nucleus removal energy (`HitNucleonRequest` with 2 nucleons) |
| 3 | `kinematics/MECKinematics`, `kinematics/MECCurrent` | Lepton and two-nucleon kinematics; event weight = fully differential 2p2h integrand with the relativistic MEC (direct − exchange, all spins) |
| tool | `apps/nuflow-mec-response` | The five CC response functions $R^{CC,CL,LL,T,T'}(q,\omega)$ and the 1p1h RFG responses, for validation |

The physics model is the relativistic MEC model of Ruiz Simo, Amaro, Barbaro, De Pace, Caballero and
Donnelly [1] (RS17 below), the weak extension of the electromagnetic MEC model of De Pace *et al.* [2]
in the relativistic Fermi gas. It is the microscopic 2p2h model that underlies the SuSAv2-MEC
parametrisations [3,4]. It was implemented from the equations in [1] and [5]; no generator code was
used, and **no parameter is fitted to the MINERvA corpus**. All couplings and form factors are those of
[1,2], and the N–Δ axial vertex is the same as in the `CCRES` channel ([pion_production.md](pion_production.md)).

---

## 1. Cross section

### 1.1 Lepton and hadron tensors

With $k=(E_\nu,\mathbf k)$, $k'$ the outgoing lepton and $Q=k-k'=(\omega,\mathbf q)$, $q=|\mathbf q|$,
the inclusive cross section per nucleus is

$$
\frac{d^2\sigma}{d\omega\,dq} = \frac{G_F^2V_{ud}^2\,q}{16\pi E_\nu^2}\;L_{\mu\nu}W^{\mu\nu},
\qquad
L^{\mu\nu}=\sum_{\rm spins} l^\mu l^{\nu*},\quad l^\mu = \bar u(k')\gamma^\mu(1-\gamma_5)u(k)
$$

(for ν̄: $\bar v(k)\gamma^\mu(1-\gamma_5)v(k')$). The full contraction keeps all ten components of the
tensors, so the vector–axial interference (the $R^{T'}$ response, sign-flipped between ν and ν̄) and
the lepton-mass terms are exact. In the frame with $\mathbf q\parallel z$ the contraction reduces to
RS17 Eq. (2) with the five responses

$$
R^{CC}=W^{00},\; R^{CL}=-\tfrac12{\rm Re}(W^{03}+W^{30}),\; R^{LL}=W^{33},\;
R^{T}=W^{11}+W^{22},\; R^{T'}=-\tfrac{i}{2}(W^{12}-W^{21}).
$$

The 2p2h hadron tensor of the RFG (RS17 Eq. 9), normalised to a nucleus with $Z$ protons
($V=3\pi^2Z/k_F^3$), is

$$
W^{\mu\nu}=\frac{V}{(2\pi)^9}\int d^3p'_1\,d^3h_1\,d^3h_2\,
\frac{M^4}{E_1E_2E'_1E'_2}\,
\delta(E'_1+E'_2-E_1-E_2-\omega')\,
\Theta(p'_1,p'_2,h_1,h_2)\; r^{\mu\nu}(\mathbf p'_1,\mathbf p'_2,\mathbf h_1,\mathbf h_2),
$$

with $\mathbf p'_2=\mathbf h_1+\mathbf h_2+\mathbf q-\mathbf p'_1$, $\Theta$ the Pauli factor
$\theta(k_F-h_1)\theta(k_F-h_2)\theta(p'_1-k_F)\theta(p'_2-k_F)$ and
$\omega'=\omega-\epsilon$ the energy transfer reduced by the pair separation energy. The elementary
tensor is $r^{\mu\nu}=2\sum_{\rm spins}A^{\mu*}A^\nu$ (RS17 Eqs. 57–63, covariant spinors
$\bar uu=2M$ absorbed into the $M^4/(E\ldots)$ factor), with the **antisymmetrised** amplitude
$A^\mu = j^\mu({\rm direct}) - j^\mu({\rm exchange})$.

### 1.2 The two-body current

The MEC is the weak pion-production amplitude of Hernández, Nieves and Valverde [5] without the
nucleon-pole terms, with the pion absorbed by a second nucleon (RS17 Eqs. 13–48):

- **seagull** (contact, vector + axial), RS17 Eqs. (17)–(18);
- **pion-in-flight** (vector), Eq. (26);
- **pion-pole** (axial, $\propto Q^\mu/(Q^2-m_\pi^2)$), Eq. (31);
- **Δ(1232) forward and backward** (vector $C_3^V$ and axial $C_5^A$), Eqs. (33)–(40), with the
  Rarita–Schwinger projector.

The pion exchanged between the lines carries $V_{\pi NN}(p',h)=\bar u(p')\gamma_5\not k\,u(h)/(k^2-m_\pi^2)$ and
$f_{\pi NN}/m_\pi = g_A/2f_\pi$ (Goldberger–Treiman). The isospin structure is written (RS17 Eq. 50) as

$$
j^\mu=\tau_+(1)J_1^\mu+\tau_+(2)J_2^\mu+(I_V)_+J_3^\mu ,
$$

which gives the effective charge-channel currents of RS17 Table I:
$np\to pp$: $J_1+J_3$ (exchange: final momenta swapped); $nn\to np$: $J_2+J_3$ (exchange: initial
momenta swapped); for ν̄, $pp\to np$: $J_1-J_3$, $np\to nn$: $J_2-J_3$.

### 1.3 Parameters

| Quantity | Value | Source |
|---|---|---|
| $g_A$, $f_\pi$ | 1.26, 93 MeV | [1] |
| $f^*_{\pi N\Delta}$ | 2.13 | [1,2] |
| $C_3^V(Q^2)=C_3^V(0)/[(1+Q^2/M_V^2)^2(1+Q^2/4M_V^2)]$ | $C_3^V(0)=2.13$, $M_V=0.84$ GeV | [1,5] |
| $C_5^A(Q^2)=C_5^A(0)/(1+Q^2/M_{A\Delta}^2)^2$ | $C_5^A(0)=1.00$, $M_{A\Delta}=0.93$ GeV | fit IV of [6]; identical to `CCRES` |
| $F_1^V(Q^2)$ (seagull, pion-in-flight) | dipole/Galster isovector | [1,5] |
| $F_\rho(k^2)=1/(1-k^2/m_\rho^2)$ (axial seagull, pion pole) | $m_\rho=775.8$ MeV | [1] |
| Strong form factors $(\Lambda^2-m_\pi^2)/(\Lambda^2-k^2)$ | $\Lambda_{\pi NN}=1.3$ GeV, $\Lambda_{\pi N\Delta}=1.15$ GeV | [2] Eqs. (16)–(19), used in [1] |
| Δ propagator | **real part** of $1/(P^2-M_\Delta^2+iM_\Delta\Gamma(P^2))$, free width $\Gamma(s)$ (P-wave, $f^*$) | [1] §V |
| $k_F$, removal energy (¹²C) | 0.221 GeV, 25 MeV per nucleon | Moniz *et al.* (`data/nuclei/nuclei.yaml`; RS17 use 0.228 GeV and $\epsilon=40$ MeV) |

Only $C_3^V$ and $C_5^A$ are kept in the Δ vertex, as in [1]. Keeping only the real part of the Δ
propagator is RS17's prescription to avoid double counting real-pion emission ($\Delta\to N\pi$, our
`CCRES` channel) with the 2p2h channel; [1] stresses that the separation is model dependent and only the
sum is meaningful.

### 1.4 Nuclear embedding and normalisation in the generator

The `NuclearModel` (RFG) gives two holes $\mathbf h_1,\mathbf h_2$, each uniform in its Fermi sphere,
with off-shell energies $E_i=\sqrt{h_i^2+M^2}-E_b$. The pair type is drawn with probability ½ each
(`SampleRequest`) and the channel multiplicity is $2N$ (ν) or $2Z$ (ν̄), so each charge channel enters
with weight $N$ (resp. $Z$), as in RS17's symmetric-matter normalisation. The per-pair density is
$P_X=k_F^3/(12\pi^2)$. With these factors the MC estimator below is identical to $d^2\sigma$ of §1.1 with
$V(4\pi k_F^3/3)^2/(2\pi)^9 = Zk_F^3/(96\pi^5)$; this was checked numerically (§3.2).

The removal energy enters as in our CCQE (de Forest-type): the kinematics conserve four-momentum with
the off-shell holes (so the pair receives $\omega - 2E_b$ of excitation), while the current is evaluated
with on-shell spinors and the effective transfer $\tilde Q = p'_1+p'_2-\tilde h_1-\tilde h_2$.

---

## 2. Sampling

`MECKinematics::Generate` samples, per trial,

1. $q$ from a proposal $\propto(q+q_s)^{-2}$ on $[E_\nu-k'_{\max},\,E_\nu+k'_{\max}]$ ($q_s=1$ GeV,
   option `proposal_q_scale`), and an optional upper cut `q3_max` (default none);
2. $\omega$ uniform between the 2p2h threshold for this pair, $\sqrt{(2M)^2+(q-|\mathbf h_1+\mathbf h_2|)^2_+}-E_{\rm pair}$, and the lepton kinematic limit;
3. the lepton azimuth uniformly;
4. the direction of $\mathbf p'_1$ isotropically in the rest frame of $P=Q+h_1+h_2$, with
   $\mathbf p'_2 = \mathbf P-\mathbf p'_1$.

Pauli blocking ($p'_{1,2}>k_F$ of their species) is applied with the same `NuclearModel::PauliBlocked`
as CCQE. The weight is the integrand

$$
w = P_X\,\frac{G_F^2V_{ud}^2\;q\,p^*\;\sum_{\rm spins}L_{\mu\nu}A^{\mu*}A^\nu}
{256\pi^3E_\nu^2\,E_1E_2\,W_{\rm pair}\;g(q,\omega)}
$$

($p^*$ and $W_{\rm pair}$ the pair CM momentum and mass, $g$ the proposal density). Because the
isotropic CM direction is only the *proposal* and the matrix element is evaluated for each
configuration, the nucleon angular and energy sharing follow the MEC model instead of the usual
isotropic-decay approximation (cf. Sobczyk [7], who discusses the limitations of that approximation).
Accepted events are unweighted by the framework (`framework.md`). The event record has the two struck
nucleons (status 21), the lepton, the two outgoing nucleons (status 22 → status 1 through stage 4) and
the remnant (status 23); four-momentum is conserved exactly. The table build (150 energy nodes ×
2×10⁵ samples) takes about 30 s on 16 threads, and generation in the loop config runs at about 1.2×10⁵ events/s
(all channels).

---

## 3. Validation

### 3.1 Response functions against RS17 (¹²C, $k_F=0.228$ GeV, $\epsilon=40$ MeV)

`nuflow-mec-response` evaluates $W^{\mu\nu}$ by MC in the RS17 conventions. Values read off RS17 Figs. 8,
10 and 11:

| $q$ [MeV] | $\omega$ [MeV] | quantity | RS17 | nuflow |
|---|---|---|---|---|
| 800 | 500 | $R^T$ / $R^{T'}$ / $R^{CC}$ [GeV⁻¹] | ≈15.3 / 7.5 / 2.6 | 16.0 / 7.7 / 2.3 |
| 800 | 640 (dip) | $R^T$ | ≈6.3 | 6.3 |
| 400 | peak (310 / 350) | $R^T$ | ≈51 | 54 |
| 400 | 300 | direct-only / (direct+exchange), $R^T$ | 71/51 = 1.39 | 1.39 (unit test) |
| 1200 | ≈790 (Δ bump) | $R^T$ | ≈3.8 | 3.9 |
| 200 | 200 (light cone) | $R^T$ | ≈33 | 21 |

The Δ-region shape (MEC peak moving as $\omega\approx\sqrt{M_\Delta^2+q^2}-M$, the dip above it, the
steep rise of $R^{CC}$ towards the light cone) is reproduced. The largest differences are at the
light-cone end-points, where the MC integration and RS17's plotted curves are both steep. The pp ($np\to pp$)
to np transverse ratio at $q=800$, $\omega=500$ is ≈6 (≈4.5 on the Δ peak at $\omega=550$), consistent with the factor "6 or less" quoted in [8].
The vector Δ current is conserved to machine precision (unit test).

### 3.2 Generator normalisation

$d^2\sigma/d\omega\,dq$ from the event weights versus the direct contraction $L_{\mu\nu}W^{\mu\nu}$ of §1.1
($E_\nu=3$ GeV, ¹²C, [10⁻³⁸ cm²/GeV² per nucleus]): 6.74±0.05 vs 6.65 ($\omega=0.5$, $q=0.8$),
5.88±0.03 vs 5.76 ($0.3$, $0.5$), 2.34±0.03 vs 2.41 ($0.8$, $1.2$). The 1–3% differences come from
the off-shell-hole versus $\omega-\epsilon$ treatment of the separation energy.

### 3.3 Total cross sections (¹²C, per nucleus, 10⁻³⁸ cm²)

| $E_\nu$ [GeV] | 0.57 | 0.79 | 1.09 | 1.51 | 2.08 | 2.88 | 3.98 | 5.49 | 10.5 | 20 |
|---|---|---|---|---|---|---|---|---|---|---|
| ν 2p2h | 0.58 | 1.07 | 1.38 | 1.53 | 1.63 | 1.65 | 1.70 | 1.75 | 1.77 | 1.81 |
| ν CCQE (same RFG) | 4.48 | 5.30 | 5.69 | 5.76 | 5.66 | 5.53 | 5.39 | 5.28 | 5.15 | 5.07 |

ν̄: 0.077 (0.5 GeV), 0.39 (1.05), 0.87 (2.19), 1.27 (4.57), 1.52 (9.56), 1.67 (20). The ν cross
section saturates above ≈2 GeV at 1.7–1.8×10⁻³⁸ cm² per ¹²C (0.15×10⁻³⁸ per nucleon, ≈33% of CCQE),
with a slow residual rise from the large-$q$ light-cone region; ν̄ rises more slowly (destructive
$R^{T'}$ interference) and approaches ν at high energy.

**MINERvA ME FHC flux on CH** (loop config, 2M events): $\langle\sigma_{\rm 2p2h}\rangle = 1.32\times10^{-39}$ cm² per
nucleon (CCQE 4.07×10⁻³⁹, ratio 0.32); 3.5% of all events. Final states: pp 77%, np 23% (the $np$
initial pair dominates, as in [1,8]). $q_3$ distribution of the 2p2h events: 8% below 0.4 GeV, 39% in
0.4–0.8, 28% in 0.8–1.2, 18% in 1.2–2, 8% above 2 GeV.

---

## 4. Limitations

- **High-$q$ extrapolation.** RS17 present responses up to $q=2$ GeV; the model is relativistic and
  we apply **no $q$ cutoff** (Valencia-model applications use $q<1.2$ GeV [9]). At the ME flux 26% of the
  2p2h cross section has $q_3>1.2$ GeV and 8% has $q_3>2$ GeV. These events lie outside the published
  validation range; the option `q3_max` exists but is off.
- **Δ-propagator prescription.** Keeping only the real part of the free Δ propagator is RS17's
  choice; it puts the MEC peak at the Δ position in $\omega$ and leaves sizable strength near the
  light cone at $q\lesssim0.4$ GeV, where the Δ pole is reached through Fermi motion. No in-medium Δ
  self-energy (Oset–Salcedo) is included. Valencia-type models with a dressed Δ put more strength in
  the QE–Δ dip [9]; the 2p2h/Δ split is model dependent [1].
- **Only MEC.** NN-correlation (one-body current with correlated pairs) and correlation–MEC
  interference terms of Valencia-type models are not included; SuSAv2 also has only MEC.
- **Nuclear model.** Global RFG ($k_F=0.221$ GeV, $E_b=25$ MeV per hole) for holes and Pauli blocking;
  no LFG, no spectral function, no short-range-correlation pair momentum.
- **Isospin fraction.** The $np$:$nn$ mixture comes from the model's matrix elements (symmetric
  matter $N=Z$); for $N\neq Z$ nuclei the channel weights use $N$ (ν) or $Z$ (ν̄) with the species
  Fermi momentum, a simple extension of RS17.
- **No FSI** of the outgoing nucleons ([fsi_none.md](fsi_none.md)).
- **Only $\nu_\mu$/$\bar\nu_\mu$ CC**; NC 2p2h is not implemented.

## 5. Configuration

```yaml
  - type: CCMEC
    nuclear_model: {type: RFG}
    mec: {C5A0: 1.00, MA: 0.93, C3V0: 2.13, MV: 0.84, f_star: 2.13, delta_propagator: real,   # or full
          lambda_piNN: 1.3, lambda_piNDelta: 1.15, mass: 1.232, gA: 1.26, f_pi: 0.093}
    nucleon_form_factors: {vector: dipole, axial: dipole}   # F1V only
    proposal_q_scale: 1.0     # GeV
    q3_max: 0                 # GeV, 0 = no cut
```

## References

1. I. Ruiz Simo, J.E. Amaro, M.B. Barbaro, A. De Pace, J.A. Caballero, T.W. Donnelly, "Relativistic model of 2p-2h meson exchange currents in (anti)neutrino scattering", *J. Phys. G* **44**, 065105 (2017) [arXiv:1604.08423].
2. A. De Pace, M. Nardi, W.M. Alberico, T.W. Donnelly, A. Molinari, "The 2p-2h electromagnetic response in the quasielastic peak and beyond", *Nucl. Phys. A* **726**, 303 (2003) [nucl-th/0304084].
3. G.D. Megias *et al.*, "Meson-exchange currents and quasielastic neutrino cross sections in the SuperScaling Approximation model", *Phys. Rev. D* **91**, 073004 (2015) [arXiv:1412.1822].
4. G.D. Megias *et al.*, "Inclusive electron scattering within the SuSAv2 meson-exchange current approach" / "Charged-current neutrino-nucleus reactions within the SuSAv2-MEC approach", *Phys. Rev. D* **94**, 093004 (2016) [arXiv:1607.08565].
5. E. Hernández, J. Nieves, M. Valverde, "Weak pion production off the nucleon", *Phys. Rev. D* **76**, 033005 (2007) [hep-ph/0701149].
6. E. Hernández, J. Nieves, M. Valverde, M.J. Vicente Vacas, *Phys. Rev. D* **81**, 085046 (2010) [arXiv:1001.4416].
7. J.T. Sobczyk, "Multinucleon ejection model for two body current neutrino interactions", *Phys. Rev. C* **86**, 015504 (2012) [arXiv:1201.3673].
8. I. Ruiz Simo *et al.*, "Emission of neutron–proton and proton–proton pairs in neutrino scattering", *Phys. Lett. B* **762**, 124 (2016) [arXiv:1607.08451].
9. J. Nieves, I. Ruiz Simo, M.J. Vicente Vacas, *Phys. Rev. C* **83**, 045501 (2011) [arXiv:1102.2777]; R. Gran *et al.*, *Phys. Rev. D* **88**, 113007 (2013) [arXiv:1307.8105].
