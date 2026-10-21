# Stage 4: `TransparentNucleus` (no final-state interactions)

`fsi/TransparentNucleus` is the trivial stage-4 model: the nucleus is treated as transparent to the hadrons made at the primary vertex.

**Algorithm.** For every particle with status 22 ("hadron before FSI"), a copy with status 1 ("final state") is appended. Its `mother` index points to the status-22 particle. The spectator remnant (status 23) is left as it is.

**Physics content.** None. The final-state hadron kinematics are the primary-vertex kinematics:
- no rescattering, absorption or charge exchange;
- no pion production;
- no nuclear optical potential;
- no knock-out of extra nucleons.

The struck nucleon also leaves without paying any further energy cost. The removal energy $\epsilon_b$ of the Fermi gas (see [quasielastic.md](quasielastic.md)) is accounted for in the remnant's excitation energy, so four-momentum is conserved exactly event by event.

**Where it matters.** Real nuclei change the observed hadronic final states considerably. In carbon at GeV energies, for example, an outgoing proton reinteracts with probability of order 30–40%. Pions from the CCRES channel ([pion_production.md](pion_production.md)) and the hadrons from CCDIS ([dis.md](dis.md)) also escape unchanged: there is no pion absorption (which would move CC1π events into the QE-like topology), charge exchange, or rescattering (which softens the $T_\pi$ spectrum), so the predicted CC1π$^+$ $T_\pi$ spectrum is harder than measured. Comparisons with data on hadron kinematics should therefore treat this model as a baseline (the "no-FSI" truth used for unfolding studies). The lepton kinematics are unaffected by FSI in the impulse approximation.

**Interface for future models.** An intranuclear cascade (INC) model derives from `FSIModel` and implements `Propagate(nucleus, event, rng)`.
- It reads the status-22 hadrons and their positions. Positions are not yet in the event record; a cascade will add a vertex position to `Particle`.
- It transports them through the nuclear density of the remnant.
- It appends the escaping particles with status 1, plus any intermediate particles with generator-specific statuses (≥ 24), declared in the NuHepMC header.

**Superseded in the loop configuration** by the intranuclear cascade [fsi_cascade.md](fsi_cascade.md) (`fsi: {type: cascade}`), which implements the interface above: `Particle::x` positions, status-1 escaping hadrons, status-25 post-FSI remnant. `TransparentNucleus` remains available as `fsi: {type: none}` (no-FSI truth).
