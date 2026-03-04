using Oceananigans
using CairoMakie
using Printf

casename = "wave_averaged"
filename = "near_inertial_$(casename)_slices.jld2"

w_xz_ts = FieldTimeSeries(filename, "w_xz")
w_xy_ts = FieldTimeSeries(filename, "w_xy")

times = w_xz_ts.times
Nt = length(times)

f = 1e-4
T_inertial = 2π / f

# Color limits from a mid-run snapshot
n_mid = Nt ÷ 2
wmax = maximum(abs, interior(w_xz_ts[n_mid])) * 0.8
wlims = (-wmax, wmax)

x_xz = xnodes(w_xz_ts)
z_xz = znodes(w_xz_ts)
x_xy = xnodes(w_xy_ts)
y_xy = ynodes(w_xy_ts)

fig = Figure(size = (1200, 800))

time_label = Observable("")
Label(fig[0, 1], time_label, fontsize = 20, tellwidth = false)

ax_xz = Axis(fig[1, 1], xlabel = "x (m)", ylabel = "z (m)", title = "x-z slice (y = Ly/2)")
ax_xy = Axis(fig[2, 1], xlabel = "x (m)", ylabel = "y (m)", title = "x-y slice (near surface)")

xz_data = Observable(interior(w_xz_ts[1], :, 1, :))
xy_data = Observable(interior(w_xy_ts[1], :, :, 1))

hm1 = heatmap!(ax_xz, x_xz, z_xz, xz_data, colormap = :balance, colorrange = wlims)
hm2 = heatmap!(ax_xy, x_xy, y_xy, xy_data, colormap = :balance, colorrange = wlims)

Colorbar(fig[3, 1], hm1, label = "w (m/s)", vertical = false, flipaxis = false)

@info "Animating $Nt frames..."

record(fig, "vertical_velocity_$(casename).mp4", 1:Nt; framerate = 12) do n
    t = times[n]
    xz_data[] = interior(w_xz_ts[n], :, 1, :)
    xy_data[] = interior(w_xy_ts[n], :, :, 1)
    time_label[] = @sprintf("w — %s   t = %.2f inertial periods", casename, t / T_inertial)
end

@info "Done! Saved vertical_velocity_$(casename).mp4"
