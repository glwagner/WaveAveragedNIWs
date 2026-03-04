using JLD2, CairoMakie, Statistics

cases = [
    ("wave_averaged",       "uᴸ (wave-averaged)",            3, 0.5),
    ("no_waves_lagrangian", "uᴬ (Lagrangian-mean hypothesis)", 1.5, 1.0),
    ("no_waves_eulerian",   "uᴬ (Eulerian-mean hypothesis)",   1.5, 1.0),
]

f = 1e-4
T_inertial = 2π / f

fig = Figure(size = (1000, 700))

ax_us = Axis(fig[1, 1], ylabel = "U (m/s)", title = "Near-surface", xticklabelsvisible = false)
ax_vs = Axis(fig[2, 1], ylabel = "V (m/s)", xlabel = "Time (inertial periods)")
ax_um = Axis(fig[1, 2], title = "Domain average", xticklabelsvisible = false)
ax_vm = Axis(fig[2, 2], xlabel = "Time (inertial periods)")

for (casename, label, lw, alpha) in cases
    file = jldopen("near_inertial_$(casename)_averages.jld2")
    iters = parse.(Int, keys(file["timeseries/t"]))
    sort!(iters)

    t = [file["timeseries/t/$i"] for i in iters]
    U_surface = [file["timeseries/U/$i"][1, 1, end-5] for i in iters]
    V_surface = [file["timeseries/V/$i"][1, 1, end-5] for i in iters]
    U_mean = [mean(file["timeseries/U/$i"]) for i in iters]
    V_mean = [mean(file["timeseries/V/$i"]) for i in iters]
    close(file)

    t_ip = t ./ T_inertial

    kw = (; label, linewidth = lw, alpha)
    lines!(ax_us, t_ip, U_surface; kw...)
    lines!(ax_vs, t_ip, V_surface; kw...)
    lines!(ax_um, t_ip, U_mean; kw...)
    lines!(ax_vm, t_ip, V_mean; kw...)
end

axislegend(ax_us, position = :rt)
save("surface_velocity_timeseries.png", fig, px_per_unit = 2)
