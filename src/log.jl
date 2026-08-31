# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

# This file has logging code.

import Dates
import ProgressMeter
import Printf

_disable_progress_bar = false

function _print_simulation_started()
    @info "Simulation started."
end

function _print_start_time()
    time = Dates.now()
    @info "Solver started at: " * Dates.format(time, "HH\\h MM\\m SS\\s")
end

function _print_finish_time()
    time = Dates.now()
    @info "Solver finished at: " * Dates.format(time, "HH\\h MM\\m SS\\s")
end

mutable struct _BarState
    prog::ProgressMeter.Progress
    start_time::Float64
    total::Int
end

const _bar = Ref{Union{_BarState,Nothing}}(nothing)

function _hms(seconds::Real)
    s = max(0, round(Int, seconds))
    h = div(s, 3600)
    m = div(mod(s, 3600), 60)
    sec = mod(s, 60)
    return Printf.@sprintf("%02dh:%02dm:%02ds", h, m, sec)
end

function _render()
    _BAR_LEN::Int = 37  # width of the "=====" fill section

    st = _bar[]
    st === nothing && return

    p = st.prog
    n = st.total
    counter = min(p.counter, n)

    pct = n == 0 ? 100 : round(Int, 100 * counter / n)
    filled = n == 0 ? _BAR_LEN : round(Int, _BAR_LEN * counter / n)
    bar = "|" * repeat("=", filled) * repeat(" ", _BAR_LEN - filled) * "|"

    elapsed = time() - st.start_time
    eta = (counter > 0 && counter < n) ? elapsed / counter * (n - counter) : 0.0

    line = Printf.@sprintf(
        "Running %s %3d%% [%s - %s]", bar, pct, _hms(elapsed), _hms(eta)
    )
    print(stdout, "\r", line)
    flush(stdout)
end

function _create_progress_bar(total_ticks::Integer)
    prog = ProgressMeter.Progress(total_ticks; enabled=false)
    _bar[] = _BarState(prog, time(), total_ticks)
    _render()
    return nothing
end

function _update_progress_bar()
    st = _bar[]
    ProgressMeter.next!(st.prog)
    _render()
    return nothing
end

function _finish_progress_bar()
    st = _bar[]
    ProgressMeter.finish!(st.prog)
    _render()
    println(stdout)
    _bar[] = nothing
    return nothing
end