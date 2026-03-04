using Oceananigans
using CairoMakie
using Printf

f = 1e-4
N² = 1e-5
h = 100
g = 9.81

# Wave parameters (MSM97)
amplitude = 0.8
wavelength = 60.0
wavenumber = 2π / wavelength
frequency = sqrt(g * wavenumber)
vertical_scale = wavelength / (4π)
Uˢ = amplitude^2 * wavenumber * frequency  # surface Stokes drift

# (Ri, Nx, Ny, Nz)
runs = [(10, 128, 128, 128),
        (1, 256, 256, 256),
        (0.1, 256, 256, 256),
        (0.01, 256, 256, 256)]

nrows = length(runs)
# Columns: b ws, b ns, uᴱ ws, uᴸ ws, u ns, w ws, w ns
fig = Figure(size = (3000, 400 * nrows + 80), fontsize = 14)

for (ri, (Ri, Nx, Ny, Nz)) in enumerate(runs)
    Ri_str = Ri >= 1 ? @sprintf("%d", Int(Ri)) : @sprintf("%g", Ri)

    # Load with_stokes data
    prefix_ws = "baroclinic_with_stokes_Ri$(Ri_str)_$(Nx)x$(Ny)x$(Nz)"
    bt_ws = FieldTimeSeries("$(prefix_ws)_top_slice.jld2", "b")
    ut_ws = FieldTimeSeries("$(prefix_ws)_top_slice.jld2", "u")
    wt_ws = FieldTimeSeries("$(prefix_ws)_top_slice.jld2", "w")

    # Load no_stokes data
    prefix_ns = "baroclinic_no_stokes_Ri$(Ri_str)_$(Nx)x$(Ny)x$(Nz)"
    bt_ns = FieldTimeSeries("$(prefix_ns)_top_slice.jld2", "b")
    ut_ns = FieldTimeSeries("$(prefix_ns)_top_slice.jld2", "u")
    wt_ns = FieldTimeSeries("$(prefix_ns)_top_slice.jld2", "w")

    # Use last snapshot
    n_ws = length(bt_ws.times)
    n_ns = length(bt_ns.times)
    t_days = bt_ws.times[n_ws] / 86400

    b_ws = interior(bt_ws[n_ws], :, :, 1)
    uL_ws = interior(ut_ws[n_ws], :, :, 1)  # uᴸ (model velocity in CL)
    uE_ws = uL_ws .- Uˢ                      # uᴱ = uᴸ - uˢ
    w_ws = interior(wt_ws[n_ws], :, :, 1)

    b_ns = interior(bt_ns[n_ns], :, :, 1)
    u_ns = interior(ut_ns[n_ns], :, :, 1)
    w_ns = interior(wt_ns[n_ns], :, :, 1)

    grid = bt_ws.grid
    x = xnodes(grid, Center()) ./ 1e3
    y = ynodes(grid, Center()) ./ 1e3

    # Data-driven u color range
    u_max = max(maximum(abs, uE_ws), maximum(abs, uL_ws), maximum(abs, u_ns), 1e-6)
    u_lim = 0.8 * u_max

    # Column 1: b with_stokes
    ax1 = Axis(fig[ri, 1],
               xlabel = ri == nrows ? "x (km)" : "",
               ylabel = @sprintf("Ri = %s\nt = %.1f d\ny (km)", Ri_str, t_days),
               title = ri == 1 ? "b (with Stokes)" : "",
               aspect = 1,
               xticklabelsvisible = ri == nrows,
               yticklabelsvisible = true)
    heatmap!(ax1, x, y, b_ws, colormap = :thermal)

    # Column 2: b no_stokes
    ax2 = Axis(fig[ri, 2],
               xlabel = ri == nrows ? "x (km)" : "",
               title = ri == 1 ? "b (no Stokes)" : "",
               aspect = 1,
               xticklabelsvisible = ri == nrows,
               yticklabelsvisible = false)
    heatmap!(ax2, x, y, b_ns, colormap = :thermal)

    # Column 3: uᴱ = uᴸ - uˢ (with_stokes)
    ax3 = Axis(fig[ri, 3],
               xlabel = ri == nrows ? "x (km)" : "",
               title = ri == 1 ? "uᴱ = uᴸ - uˢ (with Stokes)" : "",
               aspect = 1,
               xticklabelsvisible = ri == nrows,
               yticklabelsvisible = false)
    heatmap!(ax3, x, y, uE_ws, colormap = :balance, colorrange = (-u_lim, u_lim))

    # Column 4: uᴸ (with_stokes model velocity)
    ax4 = Axis(fig[ri, 4],
               xlabel = ri == nrows ? "x (km)" : "",
               title = ri == 1 ? "uᴸ (with Stokes)" : "",
               aspect = 1,
               xticklabelsvisible = ri == nrows,
               yticklabelsvisible = false)
    heatmap!(ax4, x, y, uL_ws, colormap = :balance, colorrange = (-u_lim, u_lim))

    # Column 5: u (no_stokes / wave-agnostic)
    ax5 = Axis(fig[ri, 5],
               xlabel = ri == nrows ? "x (km)" : "",
               title = ri == 1 ? "uᴬ (no Stokes)" : "",
               aspect = 1,
               xticklabelsvisible = ri == nrows,
               yticklabelsvisible = false)
    heatmap!(ax5, x, y, u_ns, colormap = :balance, colorrange = (-u_lim, u_lim))

    # Determine w color range from data
    w_max = max(maximum(abs, w_ws), maximum(abs, w_ns), 1e-6)
    w_lim = 0.8 * w_max

    # Column 6: w (with_stokes)
    ax6 = Axis(fig[ri, 6],
               xlabel = ri == nrows ? "x (km)" : "",
               title = ri == 1 ? "w (with Stokes)" : "",
               aspect = 1,
               xticklabelsvisible = ri == nrows,
               yticklabelsvisible = false)
    heatmap!(ax6, x, y, w_ws, colormap = :balance, colorrange = (-w_lim, w_lim))

    # Column 7: w (no_stokes)
    ax7 = Axis(fig[ri, 7],
               xlabel = ri == nrows ? "x (km)" : "",
               title = ri == 1 ? "w (no Stokes)" : "",
               aspect = 1,
               xticklabelsvisible = ri == nrows,
               yticklabelsvisible = false)
    heatmap!(ax7, x, y, w_ns, colormap = :balance, colorrange = (-w_lim, w_lim))
end

Label(fig[0, :], "Eulerian-mean vs Lagrangian-mean surface velocity — Uˢ = $(round(Uˢ*100, digits=1)) cm/s", fontsize = 20)

save("baroclinic_uE_vs_uL.png", fig, px_per_unit = 2)
@info "Saved baroclinic_uE_vs_uL.png"
