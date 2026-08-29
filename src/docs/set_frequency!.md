Sets the frequecies to perform the simulations with two mutually exclusive ways to use it:
- **automatic mode**: pass `max` (and optionally `min`, `length`, `sampling_density` and `step`) and Numav will automatically determine the frequency vector;
- **manual mode**: pass `vector` directly, giving full control over which frequencies are evaluated.

| Positional arguments | Type | Description |
|:--|:--|:--|
| `simulation` | `Simulation` | The simulation instance. |
| **Keyword arguments** | | |
| `max` | `Real` | Upper frequency limit (Hz). Required unless `steps` is provided. Mutually exclusive with `steps`. |
| `min` | `Real` | Lower bound of the frequency range (Hz). Only usable together with `max`. Defaults to `0` if omitted. |
| `length` | `Integer` | Number of frequency steps to compute within the `min`/`max` range. Only usable together with `max`. Defaults to `4096` if omitted. |
| `sampling_density` | `Numav.Option` | Sampling strategy to use within the range. Only usable together with `max`. Defaults to `Quadratic` if omitted. |
| `step` | `Real` | Difference in Hertz between each frequency step, supposing equally spaced steps. Only usable together with `max`. Cannot be used together with `length` and `sampling_density`. |
| `vector` | `Vector{Real}` | List of frequencies in Hertz to solve at. Cannot be used together with `max`, `min`, `length`, `sampling_density` and `step`. |

---
# `sampling_density` options

| Mode | Description |
|:--|:--|
| `Quadratic` (_default_) | Frequency discretization density grows quadratically with frequency. It tends to follow the modal density of rooms, being the recommended option for good accuracy and computation time balance. |
| `Constant` | Frequency steps are evenly spaced (uniform spacing). Suitable for broadband analyses where equal resolution at all frequencies is desired. |
 

---
# Examples

> Solve up to a maximum frequency, letting Numav automatically pick the steps from 0 Hz:
> ```julia
> set_frequency!(s, max=200) # solve from 0 Hz to 200 Hz
> ```
 
> Solve within an explicit range, by adding `min`:
> ```julia
> set_frequency!(s, min=20, max=200) # solve from 20 Hz to 200 Hz
> ```
 
> Control the number of frequency steps with `length` (defaults to `4096` when not specified):
> ```julia
> set_frequency!(s, max=200, length=500) # evaluate at 500 frequency points
> ```
> A higher step count gives finer frequency resolution at the cost of longer computation time.

> Control how steps are distributed with `sampling_density` (defaults to `Quadratic` when not specified):
> ```julia
> set_frequency!(s, max=200, sampling_density=Constant)
> ```

> Manually specify the exact frequencies with `vector`, for example to match measurement points or concentrate resolution around a resonance:
> ```julia
> set_frequency!(s, vector=[10, 40, 60, 100]) # Solve only at these frequencies
> ```

!!! tip
    Use [`get_frequency_vector`](@ref) to visually check the defined frequency vector.
