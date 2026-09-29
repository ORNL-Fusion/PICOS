# PICOS++ kinetic ECH validation (v20)

Date: 2026-09-29

## Configuration

- Initial density and ion profiles: validated hybrid steady-state average.
- Kinetic species: D+ ions and electrons, 1D-2V guiding-center push.
- Markers per cell: 330 ions and 1,320 electrons.
- Field model: reformulated Poisson with kinetic-electron quasineutral projection.
- Boundaries: current-balanced logical electron sheath and subsonic-only Bohm ion outflow.
- Collisions: self- and electron-ion collisions enabled with conservation projections.
- ECH: 300 kW, 70 GHz, harmonic 2, spatial gate 2.6--2.8 m.
- RF particle energy and velocity caps: disabled.
- Duration: 10.3861 microseconds, with matched ECH-OFF and ECH-ON initial states.

## Result

The earlier long-run failure was a marker-resolution artifact. With 82 ions/cell,
boundary current was quantized in packets equivalent to about 16 electron markers.
With 330 electrons/cell, the small electron tail crossing the downstream ambipolar
potential was under-resolved and eventually left empty cells. Increasing both
species to the counts above removed both defects.

### Stability gates

| Metric | ECH OFF | ECH ON |
|---|---:|---:|
| Final charge-density relative L2 mismatch | 0.709% | 0.660% |
| Maximum charge-density relative L2 mismatch | 1.796% | 1.796% |
| Minimum electron density over all outputs | 2.797e19 m^-3 | 3.023e19 m^-3 |
| Maximum absolute axial electric field | 39.205 V/m | 39.351 V/m |

No downstream marker hole or field runaway occurred.

### ECH gates

| Metric | Value |
|---|---:|
| Requested RF energy | 3.11584 J |
| Absorbed RF energy | 3.11768 J |
| Absorbed/requested ratio | 1.000591 |
| Final electron energy, ON minus OFF | +2.43233 J |
| Final global mean electron energy, ON minus OFF | +0.74976 eV |
| Final resonance-region perpendicular energy, OFF | 10.0119 eV |
| Final resonance-region perpendicular energy, ON | 12.8301 eV |
| Resonance-region perpendicular-energy increase | 28.15% |

The ECH-ON EEDF also contains a localized high-energy population extending to
about 200 eV near and downstream of the resonance. The additional net electron
boundary loss in the ON case is 0.37335 J, confirming that open-boundary
transport removes part of the deposited energy.

### Nuclear Fusion density gate

- Validated hybrid steady average relative L2 error: 5.56%.
- Kinetic input profile relative L2 error: 5.62%.
- Raw kinetic initial marker density relative L2 error: 8.34%.
- Raw ECH-OFF final density relative L2 error: 8.19%.

## Qualification boundary

This validates the fully kinetic ECH response for the tested 10-microsecond
window. It does not yet demonstrate a millisecond kinetic steady state. The
common OFF/ON electron cooling is the expected open-boundary energy sink relative
to the hybrid model's prescribed electron-temperature reservoir; a long steady
stage therefore still requires the explicit background-heating/thermostat model
and a separate millisecond convergence run before turning that reservoir off for
the ECH transient.
