# Conventions

---

### Complex exponential signs

The sign convention for complex acoustical equations is $+i\omega t$. Thus, sound pressure at some point in space could be represented as $p(t) = \mathrm{Re} \{ \tilde{P} e^{i\omega t} \}$, where $\tilde{P}$ is the [complex pressure amplitude](@ref "Complex amplitudes").

---

### Direction of known normal velocity conditions

To define a known normal velocity condition using [`add_sound_source!`](@ref), the velocity is considered positive when pointing **inwards**, to the direction of the simulated domain's interior.

---

### Direction of known specific acoustic impedence

Specific acoustic impedance of some surface is defined by
```math
\tilde{Z} = \frac{\tilde{P}}{\tilde{U}}
```
where:
- $ \tilde{P} $ is the [complex pressure amplitude](@ref "Complex amplitudes");
- $ \tilde{U} $ is the [complex normal velocity amplitude](@ref "Complex amplitudes").

Given that, when defining a known specific acoustic impedance using [`add_surface_material!`](@ref), the sign convention of the normal velocity is considered positive when pointing **outwards**, to the direction away from the simulated domain.
