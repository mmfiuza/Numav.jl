# Complex amplitudes

Numav deals with physical quantities that have sinusoidal oscillations, with single harmonics. Such quantities can be represented by
```math
a(t) = A \cos(\omega t + \phi)
```
where:
- $ A $ is the amplitude;
- $ \omega $ is the angular frequency;
- $ t $ is the time;
- $ \phi $ is the phase.

In this case, the physical quantity can equivalently be represented by
```math
a(t) = \mathrm{Re} \left\{ \tilde{A} e^{i\omega t} \right\},
```
where $\tilde{A}$ is the complex amplitude, given by
```math
\tilde{A} = A e^{i\phi}.
```

This approach groups amplitude and phase into a single number, the complex amplitude.
