# Charged-current single-pion production in the Δ(1232) region (CCRES)

This note covers the `nuflow` model for weak single-pion production on bound nucleons. It has three parts:

| Stage | Class | Role |
|---|---|---|
| 1 | `channel/CCRESChannel` | Either nucleon can be struck (proton with probability $Z/A$); multiplicity $A$; process id 400 |
| 2 | `initial_state/RelativisticFermiGas`, `FreeNucleon` | Struck-nucleon momentum and removal energy (same as CCQE) |
| 3 | `kinematics/ResonantPionKinematics`, `kinematics/PionProductionCurrent`, `core/Dirac.h` | Three-body final state ($\ell$, $N'$, $\pi$) and the cross-section weight |

The reactions are
$$
\nu_\ell\,p\to\ell^-p\pi^+,\qquad \nu_\ell\,n\to\ell^-n\pi^+,\qquad \nu_\ell\,n\to\ell^-p\pi^0,
$$
and for antineutrinos $\bar\nu_\ell\,n\to\ell^+n\pi^-$, $\bar\nu_\ell\,p\to\ell^+p\pi^-$, $\bar\nu_\ell\,p\to\ell^+n\pi^0$.

The hadronic current is the model of Hernández, Nieves and Valverde (HNV) [1]. It has the direct and crossed Δ(1232) pole terms and the non-resonant background required by chiral symmetry at leading order. The N–Δ axial coupling comes from the refit of Hernández *et al.* [2] to the ANL and BNL deuterium data. The amplitude is evaluated numerically with explicit Dirac matrices, and the full $\ell N\pi$ phase space is sampled. The pion angular distribution, including Δ–background interference, therefore comes from the matrix element: no decay is assumed isotropic. Everything here was implemented from the papers, not from generator code.

---

## 1. Cross section

The interaction is the contact interaction of [quasielastic.md](quasielastic.md):
$$
\mathcal M = \frac{G_F V_{ud}}{\sqrt2}\,\ell_\mu\,J^\mu,\qquad J^\mu = \bar u(p')\,S^\mu\,u(p),
$$
where $S^\mu$ is a Dirac matrix for each Lorentz index (Eq. 51 of [1], with $\cos\theta_C$ taken out). Summed over final spins and averaged over the initial nucleon spin,
$$
\overline{|\mathcal M|^2} = \frac{G_F^2V_{ud}^2}{4}\,L_{\mu\nu}H^{\mu\nu},
$$
$$
L^{\mu\nu}=\mathrm{Tr}\!\left[(\not k'+m_\ell)\gamma^\mu(1-\gamma_5)\not k\,\gamma^\nu(1-\gamma_5)\right],\qquad
H^{\mu\nu}=\mathrm{Tr}\!\left[(\not p'+M')\,S^\mu\,(\not p+M)\,\gamma^0S^{\nu\dagger}\gamma^0\right].
$$
For $\bar\nu$ the lepton line is $\mathrm{Tr}[\not k\gamma^\mu(1-\gamma_5)(\not k'-m_\ell)\gamma^\nu(1-\gamma_5)]$. Both traces are computed numerically (`PionProductionCurrent::Contract`) with the gamma matrices of `core/Dirac.h`.

The overall sign of the axial current is that of the CCQE model ($V-A$ with $g_A>0$). A test checks it: feeding the CCQE nucleon current into `Contract` reproduces the validated closed-form CCQE contraction for both $\nu$ and $\bar\nu$ (§5).

For a nucleon moving in the Fermi gas, the flux factor is $4E_\nu E_N$, the same as for CCQE (on-shell $E_N$):
$$
d\sigma=\frac{\overline{|\mathcal M|^2}}{4E_\nu E_N}\,d\Phi_3,\qquad
d\Phi_3=\frac{|\mathbf p^*_\ell|}{16\pi^2\sqrt s}d\Omega^*_\ell\;\frac{dW^2}{2\pi}\;\frac{|\mathbf p^*_\pi|}{16\pi^2W}d\Omega^*_\pi .
\tag{1}
$$
Here $s=(k+p_\text{bound})^2$, and $W$ is the $\pi N'$ invariant mass. The starred quantities are in the $k+p$ centre-of-mass frame (lepton) and in the $\pi N'$ rest frame (pion).

## 2. The HNV current (Eq. 51 of [1])

There are seven amplitudes: direct Δ pole (ΔP), crossed Δ pole (CΔP), nucleon pole (NP), crossed nucleon pole (CNP), contact term (CT), pion pole (PP) and pion-in-flight (PF). Take $q=p'+k_\pi-p$ and $f_\pi=92.4$ MeV. Then
$$
\begin{aligned}
S^\mu_{\Delta P}&= iC^{\Delta}\frac{f^*}{m_\pi}\sqrt3\,\frac{k_\pi^\alpha P_{\alpha\beta}(p+q)\,\Gamma^{\beta\mu}(p,q)}{(p+q)^2-M_\Delta^2+iM_\Delta\Gamma_\Delta}, &
S^\mu_{C\Delta P}&= iC^{C\Delta}\frac{f^*}{m_\pi\sqrt3}\,\frac{\hat\Gamma^{\mu\alpha}(p',q)P_{\alpha\beta}(p'-q)\,k_\pi^\beta}{(p'-q)^2-M_\Delta^2+iM_\Delta\Gamma_\Delta},\\
S^\mu_{NP}&=-iC^{NP}\frac{g_A}{\sqrt2f_\pi}\not k_\pi\gamma_5\frac{\not p+\not q+M}{(p+q)^2-M^2}\left[V^\mu_N-A^\mu_N\right], &
S^\mu_{CNP}&=-iC^{CNP}\frac{g_A}{\sqrt2f_\pi}\left[V^\mu_N-A^\mu_N\right]\frac{\not p'-\not q+M}{(p'-q)^2-M^2}\not k_\pi\gamma_5,\\
S^\mu_{CT}&=-iC^{CT}\frac{1}{\sqrt2f_\pi}\gamma^\mu\left(g_AF^V_{CT}\gamma_5-F_\rho(t)\right), &
S^\mu_{PP}&=-iC^{PP}F_\rho(t)\frac{1}{\sqrt2f_\pi}\frac{q^\mu}{q^2-m_\pi^2}\not q,\\
S^\mu_{PF}&=-iC^{PF}F_{PF}\frac{g_A}{\sqrt2f_\pi}\frac{(2k_\pi-q)^\mu}{t-m_\pi^2}\,2M\gamma_5, &
t&=(q-k_\pi)^2 .
\end{aligned}
$$
In these terms:
- $\hat\Gamma^{\mu\alpha}(p',q)=\gamma^0[\Gamma^{\alpha\mu}(p',-q)]^\dagger\gamma^0$.
- $F_\rho(t)=1/(1-t/m_\rho^2)$ with $m_\rho=0.7758$ GeV.
- $F_{PF}=F^V_{CT}=F_1^p-F_1^n$, which is required by CVC.
- The nucleon currents (Eq. 34 of [1]) are
$$V_N^\mu=F_1\gamma^\mu+\frac{iF_2}{2M}\sigma^{\mu\nu}q_\nu,\qquad A_N^\mu=G_A\left(\gamma^\mu\gamma_5+\frac{q^\mu\not q\,\gamma_5}{m_\pi^2-q^2}\right),$$
  with the isovector $F_{1,2}$ of the CCQE note.

**Δ ingredients.** The spin-3/2 projector (Eq. 44) is
$$
P_{\alpha\beta}(p_\Delta)=-(\not p_\Delta+M_\Delta)\left[g_{\alpha\beta}-\tfrac13\gamma_\alpha\gamma_\beta-\tfrac{2}{3}\frac{p_{\Delta\alpha}p_{\Delta\beta}}{M_\Delta^2}+\tfrac13\frac{p_{\Delta\alpha}\gamma_\beta-p_{\Delta\beta}\gamma_\alpha}{M_\Delta}\right].
$$
The width (Eq. 45) is
$$
\Gamma_\Delta(s)=\frac{1}{6\pi}\left(\frac{f^*}{m_\pi}\right)^2\frac{M}{\sqrt s}\,|\mathbf p^*_\pi(s)|^3\,\Theta(\sqrt s-M-m_\pi).
$$
It is zero in the crossed (u-channel) term. The N–Δ vertex (Eq. 40) is
$$
\Gamma^{\alpha\mu}=\Big[\tfrac{C_3^V}{M}(g^{\alpha\mu}\not q-q^\alpha\gamma^\mu)+\tfrac{C_4^V}{M^2}(g^{\alpha\mu}q\!\cdot\!p_\Delta-q^\alpha p_\Delta^\mu)+\tfrac{C_5^V}{M^2}(g^{\alpha\mu}q\!\cdot\!p-q^\alpha p^\mu)\Big]\gamma_5
+\tfrac{C_3^A}{M}(\dots)+\tfrac{C_4^A}{M^2}(g^{\alpha\mu}q\!\cdot\!p_\Delta-q^\alpha p_\Delta^\mu)+C_5^Ag^{\alpha\mu}+\tfrac{C_6^A}{M^2}q^\alpha q^\mu .
$$

**Isospin.** The coefficients $(C^\Delta,C^{C\Delta},C^{NP},C^{CNP},C^{CT},C^{PP},C^{PF})$ are $(1,1,0,1,1,1,1)$ for $\langle p\pi^+|j|p\rangle$ and $(\tfrac13,3,1,0,-1,-1,-1)$ for $\langle n\pi^+|j|n\rangle$. The $\pi^0$ amplitude is $\langle p\pi^0|j|n\rangle=-(\langle p\pi^+|j|p\rangle-\langle n\pi^+|j|n\rangle)/\sqrt2$ (Eq. 10). The antineutrino amplitudes follow from Eqs. (11)–(13). With the Δ alone the three ν channels are in the ratio $1:\tfrac19:\tfrac29$.

**Implementation.** The Δ terms are not coded as literal 4×4×4 index sums. The contractions $k_\pi^\alpha P_{\alpha\beta}$ are done analytically, so each term needs only five Dirac matrices ($X^\mu$ and $X\!\cdot\!q$). The unit test compares this with a literal transcription of Eqs. (40), (44) and (51).

## 3. Parameters

| Quantity | Value | Source |
|---|---|---|
| $C_5^A(0)$, $M_A^\Delta$ (dipole) | 1.00, 0.93 GeV | Fit IV of [2]: ANL $d\sigma/dQ^2$ ($W<1.4$) and BNL $\sigma(E_\nu\le1.1)$, with deuteron effects, flux errors and background terms |
| $C_4^A$, $C_3^A$, $C_6^A$ | $-C_5^A/4$, 0, $C_5^AM^2/(m_\pi^2+Q^2)$ | Adler model, PCAC ([1] Eqs. 47–48) |
| $C_3^V,C_4^V,C_5^V$ at 0 | 2.13, −1.51, 0.48; $M_V=0.84$ GeV; modified dipoles of [1] Eq. (46) | Electroproduction helicity amplitudes (Lalakulich–Paschos–Piranishvili [3]) |
| $M_\Delta$, $f^*$ | 1.232 GeV, 2.14 | [1] |
| $g_A$, $M_A$ (nucleon, background) | 1.26, 1.05 GeV | [1] Eq. (37) |
| Nucleon vector form factors (background) | dipole $M_V^2=0.71$ + Galster | [1] Eq. (35) |
| $m_\pi$ inside amplitudes | isospin average 138.0 MeV | [1] footnote 3 |
| $W$ range | $M+m_\pi<W<W_\text{max}=1.6$ GeV | this work, see §6 |

For comparison, the Goldberger–Treiman value is $C_5^A(0)=\sqrt{2/3}\,f_\pi f^*/m_\pi=1.17$. None of these parameters was adjusted to the MINERvA corpus.

**Configuration** (`config/minerva_me_fhc_numu_CH_loop.yaml`):
```yaml
- type: CCRES
  nuclear_model: {type: RFG}
  delta: {C5A0: 1.00, MA: 0.93}          # also: mass, f_star
  background: true                       # false = Delta pole only
  nucleon_form_factors: {vector: dipole, axial: dipole, gA: 1.26, MA: 1.05}
  w_max: 1.6                             # GeV
  # proposal_q2_scale: 1.0, proposal_bw_fraction: 0.7  (sampling only)
```

## 4. Nuclear embedding and sampling

**Initial state.** Stage 1 draws the struck nucleon's species, a proton with probability $Z/A$, through the `InteractionChannel::SampleRequest` hook. The multiplicity is $A$, so $\sigma_A=Z\sigma_p+N\sigma_n$ exactly. On a neutron, one of the two final states ($n\pi^+$ or $p\pi^0$) is drawn with probability ½ and its weight is doubled. The RFG momentum and removal energy are the same as for CCQE. Four-momentum is conserved with the off-shell bound nucleon $p_\text{bound}=(\sqrt{\mathbf p^2+M^2}-\epsilon_b,\mathbf p)$. The hadronic current uses the on-shell $p_\text{on}$ and $\tilde q=p'+k_\pi-p_\text{on}$ (the de Forest prescription, as for CCQE). Outgoing nucleons with $|\mathbf p'|<k_F$ are Pauli blocked.

**Proposal densities.** These affect only the sampling efficiency, not the distribution:
- $W$: a mixture of 70% Breit–Wigner around $M_\Delta$ (width 0.12 GeV, truncated) and 30% flat on $[M'+m_\pi,\min(W_\text{max},\sqrt s-m_\ell)]$.
- Lepton CM angle: $(1+Q^2/1\,\text{GeV}^2)^{-2}$ in $1-\cos\theta^*$, the same construction as CCQE.
- Pion: isotropic in the $\pi N'$ rest frame.

The weight is Eq. (1) divided by these densities. The unweighting envelope is the usual tabulated maximum ×1.2. In the 200k-event check, the largest $w/w_\text{env}$ was 0.93 (no violations) and the acceptance was 12%.

**Event record.** The record holds:
- the beam (4), the target (20) and the struck nucleon (21);
- the lepton (1);
- the nucleon and pion as status 22, which `TransparentNucleus` copies to status 1;
- the remnant (23).

$\pi^0$ is left undecayed, as the measurements count it as a meson. `W` in the event tree is the true $\pi N'$ invariant mass.

## 5. Validation

Unit tests (`tests/test_pion_production.cc`) check the following:
1. The optimised ΔP and CΔP amplitudes agree with a literal transcription of [1] Eqs. (40), (44) and (51) to $10^{-10}$.
2. **CVC.** With the axial current switched off, $q_\mu\bar u'S^\mu u=0$ to $10^{-10}$ for $p\pi^+$, $n\pi^+$ and $p\pi^0$, including the background. This checks the relative signs and form factors of the NP, CNP, CT and PF terms and the Δ vertex.
3. `Contract` equals the closed-form CCQE contraction for a nucleon current, for $\nu$ and $\bar\nu$ (conventions and normalization).
4. With the Δ alone, the channel ratios are exactly $1:\tfrac19:\tfrac29$, and $\Gamma_\Delta(M_\Delta^2)=0.11$–$0.125$ GeV.
5. The free-proton σ lies within bounds set by the bubble-chamber data. Generated events on CH conserve four-momentum exactly, contain exactly one charged lepton, one nucleon and one pion, and have $W\le W_\text{max}$.

**Free-nucleon cross sections** (`nuflow-xsec`, $W<1.6$ GeV, $10^{-38}$ cm²) compared with the ANL/BNL $\nu_\mu p\to\mu^-p\pi^+$ data as re-extracted by Wilkinson *et al.* [4] (from ratios to CCQE, no $W$ cut):

| $E_\nu$ [GeV] | 0.66 | 0.84 | 1.08 | 1.38 | 1.77 | 2.26 | 2.90 | 4.76 |
|---|---|---|---|---|---|---|---|---|
| nuflow $p\pi^+$ | 0.24 | 0.41 | 0.55 | 0.65 | 0.71 | 0.74 | 0.75 | 0.77 |
| ANL / BNL [4] (approx., read from Fig. 6) | 0.28 / 0.26 | – / 0.44 (0.9) | 0.47–0.50 / 0.58 | – | 0.71 / 0.66 | – | 0.74 / 0.82 (2.6) | 0.92 / 1.0 (4.5, ±0.1–0.35) |

The agreement is good up to 3 GeV. Above that the data include higher resonances ($W>1.6$) that the model does not contain. The free-neutron sum $\sigma(n\pi^+)+\sigma(p\pi^0)$ is 0.22, 0.33 and 0.41 at 1.0, 1.6 and 2.5 GeV, which is about half of $p\pi^+$. On carbon, 35% of neutron events are $n\pi^+$ and 65% are $p\pi^0$.

**Pion angular distribution.** For $p\pi^+$ on hydrogen in the ME flux, with $1.15<W<1.30$ GeV, the Adler-frame distribution is peaked at $\cos\theta^*_\pi\approx0$ (0.77 at $|\cos|\to1$ relative to its mean), with $\langle P_2\rangle=-0.054$. This corresponds to $\rho_{33}\approx0.64$ (normalized so that $\rho_{33}+\rho_{11}=1$, isotropic = 0.5), the expected M1-dominated, non-isotropic Δ decay.

**MINERvA ME flux** (CH, `minerva_me_fhc_numu_CH_loop.yaml`): $\langle\sigma_\text{CCRES}\rangle=0.547\times10^{-38}$ cm²/nucleon, against $\langle\sigma_\text{CCQE}\rangle=0.408$, a ratio of 1.34. 81% of CCRES events have $W<1.4$ GeV. Generating 2M events takes about 20 s, and building the tables takes about 50 s on 16 threads.

## 6. Limitations

- **$W$ range.** HNV is constructed for $W\lesssim1.4$ GeV. We generate up to $W_\text{max}=1.6$ GeV so that Fermi smearing between the true $W$ and the measured $W_\text{exp}$ (nucleon at rest) does not cut into the $W_\text{exp}<1.4$ CC1π signal. The value is a choice, not fitted to the corpus. Without higher resonances the Δ tail plus background between 1.4 and 1.6 GeV is only an approximation. The $p\pi^+$ cross section at 6 GeV is 0.66 with $W_\text{max}=1.4$, 0.78 with 1.6 and 1.03 with 2.0. The `CCDIS` channel ([dis.md](dis.md)) starts at $W_\text{min}=W_\text{max}=1.6$ GeV with the same definition of $W$, so there is no gap or overlap. With `CCNSTAR` the step in $d\sigma/dW$ at 1.6 GeV is about 1.6× on CH in the ME flux; it was 3.2× with CCRES alone.
- **Higher resonances** are a separate channel: `CCNSTAR` ([nstar.md](nstar.md)) adds $N(1440)$, $N(1520)$ and $N(1535)$ incoherently, so there is no interference with the Δ and the background. The $I=3/2$ states $\Delta(1600)$–$\Delta(1700)$ and non-resonant two-pion production below $W=1.6$ GeV are still missing. Above 1.6 GeV the inclusive strength comes from `CCDIS`, with its single-pion part fitted to ANL/BNL ([dis.md](dis.md) §3.1).
- **No unitarization.** There is no Watson-theorem phase between the Δ and the background (Alvarez-Ruso *et al.* [5]). Fits with it give a larger $C_5^A(0)$.
- **No in-medium Δ modification.** The Δ self-energy and width in the nucleus (Oset–Salcedo) and RPA are not included.
- **Pion FSI** is now applied by the cascade ([fsi_cascade.md](fsi_cascade.md)). In ME events, 19.8% of primary π⁺ on carbon are absorbed and 5.2% are charge-exchanged. Δ-peak absorption is about 2× below the π–¹²C data, so the FSI loss is probably underestimated. The Δ is produced and decays at the vertex: there is no Δ propagation and no in-medium Δ in the primary production.
- **Δ in 2p2h.** The Δ-MEC of the `CCMEC` channel ([mec.md](mec.md)) uses the same N–Δ vertex ($C_5^A(0)=1.00$, $M_A^\Delta=0.93$ GeV) but only the real part of the Δ propagator, following Ruiz Simo *et al.*, to limit double counting with real pion emission here. The split between the two channels is model dependent.
- **No coherent pion production.** Coherent production (process 100) is a separate channel.

## References

1. E. Hernández, J. Nieves, M. Valverde, "Weak pion production off the nucleon", *Phys. Rev. D* **76**, 033005 (2007) [hep-ph/0701149].
2. E. Hernández, J. Nieves, M. Valverde, M.J. Vicente Vacas, "N–Δ(1232) axial form factors from weak pion production", *Phys. Rev. D* **81**, 085046 (2010) [arXiv:1001.4416].
3. O. Lalakulich, E.A. Paschos, G. Piranishvili, *Phys. Rev. D* **74**, 014009 (2006) [hep-ph/0602210].
4. C. Wilkinson *et al.*, "Reanalysis of bubble chamber measurements of muon-neutrino induced single pion production", *Phys. Rev. D* **90**, 112017 (2014) [arXiv:1411.4482].
5. L. Alvarez-Ruso, E. Hernández, J. Nieves, M.J. Vicente Vacas, "Watson's theorem and the N Δ(1232) axial transition", *Phys. Rev. D* **93**, 014016 (2016) [arXiv:1510.06266].
6. R. González-Jiménez *et al.*, "Electroweak single-pion production off the nucleon: from threshold to high invariant masses", *Phys. Rev. D* **95**, 113007 (2017) [arXiv:1612.05511]. Discusses why the low-energy model is unphysical at high $W$; a possible route for a future $W>1.4$ extension.
