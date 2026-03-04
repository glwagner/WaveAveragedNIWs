using Oceananigans
using CairoMakie
using Printf

f = 1e-4
N² = 1e-5
h = 100

Ri_values = [1, 0.1, 0.01]
cases = [:with_stokes, :no_stokes]
Nx, Ny, Nz = 256, 256, 64

# Plot at the final snapshot for each Ri
nrows = length(Ri_values)
fig = Figure(size = (1800, 400 * nrows + 50), fontsize = 14)

for (ri, Ri) in enumerate(Ri_values)
    Ri_str = Ri >= 1 ? @sprintf("%d", Int(Ri)) : @sprintf("%g", Ri)

    for (ci, case) in enumerate(cases)
        prefix = "baroclinic_$(case)_Ri$(Ri_str)_$(Nx)x$(Ny)x$(Nz)"
        filename = "$(prefix)_top_slice.jld2"

        bt = FieldTimeSeries(filename, "b")
        ut = FieldTimeSeries(filename, "u")

        times = bt.times
        n = length(times)  # last snapshot
        t_days = times[n] / 86400

        b_data = interior(bt[n], :, :, 1)
        u_data = interior(ut[n], :, :, 1)

        grid = bt.grid
        x = xnodes(grid, Center()) ./ 1e3
        y = ynodes(grid, Center()) ./ 1e3

        # Buoyancy panel
        col_b = ci
        label = ci == 1 ? @sprintf("Ri = %s\nt = %.1f d", Ri_str, t_days) : ""
        ax_b = Axis(fig[ri, col_b],
                    xlabel = ri == nrows ? "x (km)" : "",
                    ylabel = label,
                    title = ri == 1 ? "b — $(replace(string(case), "_" => " "))" : "",
                    aspect = 1,
                    xticklabelsvisible = ri == nrows,
                    yticklabelsvisible = ci == 1)

        heatmap!(ax_b, x, y, b_data, colormap = :thermal)

        # u-velocity panel
        col_u = ci + 2
        U_g = sqrt(N² / Ri) * h
        u_lim = max(1.5 * U_g, 0.01)

        ax_u = Axis(fig[ri, col_u],
                    xlabel = ri == nrows ? "x (km)" : "",
                    ylabel = "",
                    title = ri == 1 ? "u — $(replace(string(case), "_" => " "))" : "",
                    aspect = 1,
                    xticklabelsvisible = ri == nrows,
                    yticklabelsvisible = false)

        heatmap!(ax_u, x, y, u_data, colormap = :balance, colorrange = (-u_lim, u_lim))
    end
end

Label(fig[0, :], "Baroclinic adjustment — final snapshots (30 Eady timescales)", fontsize = 20)

save("baroclinic_adjustment_final.png", fig, px_per_unit = 2)
@info "Saved baroclinic_adjustment_final.png"
