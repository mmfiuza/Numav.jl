# Copyright (c) 2026 Matheus Machado Fiuza <matheusmachadofiuza@gmail.com>

import Dates

_disable_progress_bar = false

function _print_opening_fem_helmholtz()
    @info "FEM-Helmholtz simulation started."
end

function _print_start_time()
    time = Dates.now()
    @info "Solver started at: " * Dates.format(time, "HH\\h MM\\m SS\\s")
end

function _print_finish_time()
    time = Dates.now()
    @info "Solver finished at: " * Dates.format(time, "HH\\h MM\\m SS\\s")
end

using ProgressMeter
using Printf

mutable struct _BarState
    prog::Progress
    start_time::Float64
    total::Int
end

const _bar = Ref{Union{_BarState,Nothing}}(nothing)


function _hms(seconds::Real)
    s = max(0, round(Int, seconds))
    h = div(s, 3600)
    m = div(mod(s, 3600), 60)
    sec = mod(s, 60)
    return @sprintf("%02dh:%02dm:%02ds", h, m, sec)
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

    line = @sprintf(
        "Running %s %3d%% [%s - %s]", bar, pct, _hms(elapsed), _hms(eta)
    )
    print(stdout, "\r", line)
    flush(stdout)
end


function _create_progress_bar(total_ticks::Integer)
    prog = Progress(total_ticks; enabled=false)
    _bar[] = _BarState(prog, time(), total_ticks)
    _render()
    return nothing
end


function _update_progress_bar()
    st = _bar[]
    next!(st.prog)
    _render()
    return nothing
end


function _finish_progress_bar()
    st = _bar[]
    finish!(st.prog)
    _render()
    println(stdout)
    _bar[] = nothing
    return nothing
end
