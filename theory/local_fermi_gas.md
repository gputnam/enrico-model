# Local Fermi gas (`LFG`) and Valencia RPA for the CCQE vertex

| Stage | Class | Role |
|---|---|---|
| 2 | `initial_state/LocalFermiGas` | Interaction point from the nuclear density, struck nucleon in the local Fermi sphere, Q-value energy balance, local Pauli blocking |
| 3 | `kinematics/ValenciaRPA`, `QuasiElasticKinematics` (`rpa: true`) | Long-range RPA correlations of the 1p1h response |
| tool | `apps/nuflow-lfg-validate` | RPA factors, $n(k)$, $d\sigma/dQ^2$, $\sigma(E)$, 12C$(e,e')$ for RFG / LFG / SF (± RPA) |

Both are configuration options of the `CCQE` channel. The step-17 loop configuration uses the
spectral function of [spectral_function.md](spectral_function.md) with `rpa: false`; the LFG and the
RPA remain available (and tested) because the LFG is the ground state in which the RPA is derived.
Everything is implemented from Nieves, Amaro and Valverde [1] (NAV04), not from generator code.

---

## 1. Local Fermi gas

- Interaction point $\mathbf r$ drawn from the point-nucleon density $\rho(r)$ of `fsi/NuclearDensity`
  (12C: harmonic oscillator, $b=1.590$ fm, $\rho(0)=0.179$ fm⁻³), the same density the cascade uses.
- Struck nucleon uniform in the local Fermi sphere $k<k_F^N(r)=[3\pi^2\rho_N(r)]^{1/3}$.
- Pauli blocking: outgoing nucleon needs $|\mathbf p'|>k_F^{N'}(r)$ at the same $r$.
- The interaction point is passed to the cascade (`InitialState::position`).

## 2. Energy balance

Two conventions (`energy: in_medium | asymptotic`):

- **`in_medium`** (NAV04 Eq. 43): both nucleons are on shell in a common local potential and the
  transition costs the gap $\Delta(r)=Q-Q_\text{LFG}(r)$, $Q_\text{LFG}=E_F^{p}(r)-E_F^{n}(r)$, with the
  ground-state Q value $Q=M(^{12}\text{N})-M(^{12}\text{C})=16.827$ MeV (NAV04 Table I; `Qnu`/`Qnubar` in
  `data/nuclei/nuclei.yaml`). Vertex: $E'(\mathbf k+\mathbf q)=E(\mathbf k)+q^0-\Delta(r)$. The struck
  nucleon sits in a well $V(r)=T_F(r)+S_N$, i.e. $E_R=V(r)-T(k)\in[S_N,S_N+T_F(r)]$ (mean 28.6 MeV for
  12C neutrons); the outgoing nucleon leaves a well $V(r)-\Delta(r)$ and is decelerated along its
  direction, the remnant taking the recoil.
- **`asymptotic`**: $E=\sqrt{M^2+k^2}-V(r)$ at the vertex and the outgoing nucleon is asymptotic.

## 3. Validation of the LFG

$\langle k\rangle=0.158$ GeV and $\langle E_R\rangle=28.6$ MeV (12C neutrons); vertex rms radius equal to
the point rms, $\langle k^2\rangle=\frac35\langle k_F^2(r)\rangle$ and the energy gap (16.83 + 1.29 MeV) are
unit-tested. In the step-17 trials the LFG was rejected for 12C: its $n(k)$ has a large excess at
$k<0.08$ GeV from the dilute surface (the semiclassical LDA artefact for light nuclei), which puts
2× too many events at $p_n<0.05$ GeV/c in cc0pinp_stki_me_fhc_pn (χ² 355 → 1618–1934) and makes the
12C$(e,e')$ QE peak 20–47% too high and too narrow ([spectral_function.md](spectral_function.md) §3).

---

## 4. Valencia RPA

### 4.1 Factors
In nuclear matter at the local density, the RPA sum of particle-hole (ph) and Δ-hole bubbles with the
effective spin-isospin interaction renormalises the components of the 1p1h hadron tensor (NAV04
Sec. III A, Appendix A 2). With $\mathbf q\parallel\hat z$:
$$
C_N=\frac1{|1-c_0f'(\rho)U_N|^2},\qquad C_L=\frac1{|1-U V_l|^2},\qquad C_T=\frac1{|1-U V_t|^2},
\qquad U=U_N+U_\Delta ,
$$
$$
V_l=\frac{f^2}{m_\pi^2}\Big[\Big(\frac{\Lambda_\pi^2-m_\pi^2}{\Lambda_\pi^2-q^2}\Big)^2\frac{\mathbf q^2}{q^2-m_\pi^2}+g'\Big],\quad
V_t=\frac{f^2}{m_\pi^2}\Big[C_\rho\Big(\frac{\Lambda_\rho^2-m_\rho^2}{\Lambda_\rho^2-q^2}\Big)^2\frac{\mathbf q^2}{q^2-m_\rho^2}+g'\Big],
$$
$f'(\rho)=\frac{\rho}{\rho_0}f'_\text{in}+(1-\frac{\rho}{\rho_0})f'_\text{ex}$ with $\rho_0$ the central density.
Parameters (NAV04): $g'=0.63$, $f^2/4\pi=0.08$, $\Lambda_\pi=1.2$ GeV, $C_\rho=2$, $\Lambda_\rho=2.5$ GeV,
$m_\rho=0.770$ GeV, $f^*/f=2.13$, $c_0=0.380$ GeV fm³, $f'_\text{in}=0.33$, $f'_\text{ex}=0.45$.

$U_N$ is the non-relativistic Lindhard function of symmetric matter (direct + crossed, degeneracy 4),
$U_N(z,q)=\frac{M}{\pi^2|\mathbf q|}[G(a_1)-G(a_2)]$, $G(a)=ak_F+\frac12(k_F^2-a^2)\ln\frac{a+k_F}{a-k_F}$,
$a_{1,2}=Mz/|\mathbf q|\mp|\mathbf q|/2$; $U_\Delta=\frac49(f^*/f)^2\rho\,[\frac1{z-\omega_\Delta+i\Gamma_\Delta/2}-\frac1{z+\omega_\Delta}]$,
$\omega_\Delta=\sqrt{M_\Delta^2+\mathbf q^2}-M$ (static Δ).

### 4.2 Which structure gets which factor (Table 2)
The tensor is split by form-factor structure and by block (00, 0z, zz, transverse xx+yy; mixed 0T/zT
blocks vanish after azimuthal averaging and get no factor):

| structure | 00 | 0z | zz | TT |
|---|---|---|---|---|
| $F_1^2$ | $C_N$ | $C_N$ | $C_N$ | 1 |
| $F_1F_2$ | $C_N$ | $C_N$ | $C_N$ | $C_T$ |
| $F_2^2$ | 1 | 1 | 1 | $C_T$ |
| $F_A^2$ | 1 | $C_L$ | $C_L$ | $C_T$ |
| $F_AF_P$, $F_P^2$ | $C_L$ | $C_L$ | $C_L$ | – |
| $(F_1+F_2)F_A$ (V–A, antisymmetric) | – | – | – | $C_T$ |

With all factors 1 the contraction reproduces `ContractQE` to $10^{-9}$; a brute-force Dirac-trace
decomposition checks each block (unit test).

### 4.3 Kinematic arguments and regularisation
All functions are evaluated at the ph energy $\nu=(|\mathbf k+\mathbf q|^2-k^2)/2M$, clamped into the
open 1p1h support of the local Fermi sea, and at $|\mathbf q|$. Two widths keep the factors finite
(details and justification in [spectral_function.md](spectral_function.md) §5): $\Gamma_\Delta=117$ MeV in the forward
Δh propagator, and $\Gamma_{ph}=5$ MeV in $U_N$ (complex energy $\nu+i\Gamma_{ph}/2$). Without them the
undamped factors have poles (zero-sound branch at $|\mathbf q|<0.1$ GeV, static-Δ pole at
$|\mathbf q|>0.6$ GeV) that produced weights $10^6$ times the envelope.

### 4.4 Range of validity
The non-relativistic effective interaction is applied in full for $|\mathbf q|\le0.6$ GeV and faded
linearly to no correction at $|\mathbf q|=1.0$ GeV (`q_full`, `q_zero`).

## 5. RPA results

At $\rho=0.16$ fm⁻³ along $\nu=q^2/2M$: $C_T$ = 0.13, 0.16, 0.30, 0.51, 0.92 and $C_L$ = 0.17, 0.97,
1.34, 1.27, 1.10 at $|\mathbf q|$ = 0.05, 0.2, 0.3, 0.4, 0.6 GeV. For 3 GeV $\nu_\mu$ on 12C,
$d\sigma/dQ^2$(RPA)/(no RPA) = 0.63–0.65 (LFG) and 0.70–0.75 (SF) at $Q^2<0.1$ GeV², 1.0 at
$Q^2\approx0.4$, up to 1.06 at 0.5–0.6 GeV², 1 above 0.8 GeV². $\sigma_\text{CCQE}$ is reduced by 9% (LFG) and
8% (SF) at 3–10 GeV, as in NAV04 (≈10% at 1 GeV).

## 6. Validation tool

`nuflow-lfg-validate factors | momentum | q2 E | xsec | eqe E theta`.

## 7. Configuration reference

```yaml
  - type: CCQE
    nuclear_model: {type: LFG, energy: in_medium}   # or asymptotic; or {type: SF} (spectral_function.md)
    form_factors: {vector: BBBA05, axial: dipole, MA: 1.03}
    rpa: true            # or a map with any of:
    # rpa: {enabled: true, g_prime: 0.63, f2_over_4pi: 0.08, lambda_pi: 1.2, c_rho: 2.0, lambda_rho: 2.5,
    #       m_rho: 0.770, f_star_over_f: 2.13, m_delta: 1.232, gamma_delta: 0.117, gamma_ph: 0.005,
    #       c0: 0.380, f0p_in: 0.33, f0p_ex: 0.45, delta_hole: true, q_full: 0.6, q_zero: 1.0}
```

## References

1. J. Nieves, J. E. Amaro and M. Valverde, "Inclusive quasielastic charged-current neutrino-nucleus reactions", *Phys. Rev. C* **70**, 055503 (2004), [nucl-th/0408005](https://arxiv.org/abs/nucl-th/0408005); erratum *ibid.* **72**, 019902 (2005).
2. C. García-Recio, J. Nieves and E. Oset, *Nucl. Phys. A* **547**, 473 (1992) (Lindhard functions $U_N$, $U_\Delta$).
3. W. J. Huang *et al.*, M. Wang *et al.*, *Chin. Phys. C* **45**, 030002/030003 (2021) (AME2020 masses, Q values).
