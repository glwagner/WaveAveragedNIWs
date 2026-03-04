# Baroclinic adjustment with and without wave-averaged Stokes drift
#
# Based on WaveAveragedBaroclinicAdjustment (Wagner, unpublished).
# A density front in thermal wind balance undergoes baroclinic instability.
#
# Two cases per Ri:
#   1. "with_stokes": Craik-Leibovich equations with steady Stokes drift
#   2. "no_stokes": No wave forcing (standard baroclinic adjustment)
#
# Parameter study: vary Ri (Richardson number) with fixed f, waves.
# Key scaling: Ps = f / (Ri ∂z_uˢ) ~ 1/Ri
#   Ri = 100 → Ps ~ 7e-5 (QG, waves negligible)
#   Ri = 1   → Ps ~ 7e-3 (ageostrophic, waves start to matter)
#   Ri = 0.1 → Ps ~ 0.07 (symmetric instability, waves important)

using Oceananigans
using Oceananigans.Units
using Random
using Printf
using CUDA

Random.seed!(42)

# Fixed environmental parameters
f = 1e-4        # Coriolis parameter [s⁻¹]
g = 9.81

# Wave parameters (same as McWilliams, Sullivan & Moeng 1997)
amplitude = 0.8            # wave amplitude [m]
wavelength = 60.0          # wavelength [m]
wavenumber = 2π / wavelength
frequency = sqrt(g * wavenumber)  # deep water dispersion
vertical_scale = wavelength / (4π)  # 1/(2k) decay scale
Uˢ = amplitude^2 * wavenumber * frequency  # surface Stokes drift [m s⁻¹]
∂z_uˢ_surface = Uˢ / vertical_scale

# Stokes drift
stokes_drift_parameters = (; Uˢ, h_s = vertical_scale)
@inline ∂z_uˢ(z, t, p) = p.Uˢ / p.h_s * exp(z / p.h_s)

# IVP parameters (fixed)
N² = 1e-5       # buoyancy frequency squared [s⁻²]
h = 100         # mixed layer / front depth [m]
H = 200         # total depth [m]

# Deformation radius (fixed for all Ri since N², h, f are fixed)
Ld = sqrt(N²) * h / f

# Grid resolution
Nx = 128
Ny = 128
Nz = 128

# Initial conditions
step(y, Δy) = min(max(0, y/Δy + 1/2), 1)
ramp(y, Δy=1) = max(0, y / Δy)

function run_simulation(; case, Ri, stop_time = 30days)

    # Front parameters (depend on Ri)
    M² = sqrt(f^2 * N² / Ri)
    Δy = Ld                    # front width = deformation radius
    Δb = M² * Δy              # buoyancy jump across front
    ϵb = 1e-2 * Δb

    # Derived scales
    U_g = M² * h / f          # geostrophic velocity
    Ro = U_g / (f * Ld)       # Rossby number = 1/√Ri
    Ps = f / (Ri * ∂z_uˢ_surface)

    # Domain (20 deformation radii)
    Lx = 20Ld
    Ly = 20Ld
    Lz = H

    bᵢ(x, y, z) = N² * z + Δb * step(y, Δy) * ramp(1 + z/h) + ϵb * randn()

    grid = RectilinearGrid(GPU(),
                           size = (Nx, Ny, Nz),
                           x = (0, Lx),
                           y = (-Ly/2, Ly/2),
                           z = (-Lz, 0),
                           halo = (5, 5, 5),
                           topology = (Periodic, Bounded, Bounded))

    coriolis = FPlane(; f)

    if case == :with_stokes
        stokes_drift = UniformStokesDrift(; ∂z_uˢ, parameters = stokes_drift_parameters)
    elseif case == :no_stokes
        stokes_drift = nothing
    else
        error("Unknown case: $case. Use :with_stokes or :no_stokes.")
    end

    model = NonhydrostaticModel(grid; coriolis, stokes_drift,
                                buoyancy = BuoyancyTracer(),
                                tracers = :b,
                                advection = WENO(order=9))

    Random.seed!(42)
    set!(model, b = bᵢ)

    simulation = Simulation(model; Δt = 10seconds, stop_time)
    conjure_time_step_wizard!(simulation, IterationInterval(10), cfl = 0.7, max_Δt = 10minutes)

    # Progress callback
    wall_clock = Ref(time_ns())

    function progress(sim)
        u, v, w = sim.model.velocities
        elapsed = (time_ns() - wall_clock[]) / 1e9
        @printf("Ri=%g %s | [%05.2f%%] i: %d, t: %s, wall: %s, max|u|: (%.1e, %.1e, %.1e), Δt: %s\n",
                Ri, case,
                100 * time(sim) / sim.stop_time,
                iteration(sim),
                prettytime(time(sim)),
                prettytime(elapsed),
                maximum(abs, interior(u)),
                maximum(abs, interior(v)),
                maximum(abs, interior(w)),
                prettytime(sim.Δt))
        wall_clock[] = time_ns()
    end

    add_callback!(simulation, progress, IterationInterval(100))

    # Output filenames encode Ri and resolution
    Ri_str = Ri >= 1 ? @sprintf("%d", Int(Ri)) : @sprintf("%g", Ri)
    prefix = "baroclinic_$(case)_Ri$(Ri_str)_$(Nx)x$(Ny)x$(Nz)"
    b = model.tracers.b
    u, v, w = model.velocities

    # Slices (top and east faces)
    simulation.output_writers[:top] =
        JLD2Writer(model, (; b, u, v, w),
                   schedule = TimeInterval(1day),
                   filename = "$(prefix)_top_slice.jld2",
                   indices = (:, :, Nz),
                   overwrite_existing = true)

    simulation.output_writers[:east] =
        JLD2Writer(model, (; b, u, v, w),
                   schedule = TimeInterval(1day),
                   filename = "$(prefix)_east_slice.jld2",
                   indices = (Nx, :, :),
                   overwrite_existing = true)

    # Zonal averages
    B = Average(b, dims = 1)
    U = Average(u, dims = 1)
    V = Average(v, dims = 1)
    W² = Average(w^2, dims = 1)

    simulation.output_writers[:zonal] =
        JLD2Writer(model, (; B, U, V, W²),
                   schedule = TimeInterval(1day),
                   filename = "$(prefix)_zonal_average.jld2",
                   overwrite_existing = true)

    @info "Running baroclinic adjustment: $case, Ri=$Ri for $(prettytime(stop_time))..."
    @info "  Ro = $(round(Ro, sigdigits=3)), Ps = $(round(Ps, sigdigits=3))"
    @info "  M² = $(round(M², sigdigits=3)) s⁻², U_g = $(round(U_g, sigdigits=3)) m/s"
    @info "  Ld = $(round(Ld, digits=1)) m, domain = $(round(Lx, digits=0)) × $(round(Ly, digits=0)) × $Lz m"
    @info "  Δb = $(round(Δb, sigdigits=3)) m/s², Uˢ = $(round(Uˢ, sigdigits=3)) m/s"
    run!(simulation)

    @info "Case $case Ri=$Ri complete."
    return simulation
end

# Eady timescale: τ = √Ri / (0.31 f)
eady_timescale(Ri) = sqrt(Ri) / (0.31 * f)

# (Ri, stop_time in days)
runs = [(10, 60days)]
#runs = [(1, 37days), (0.1, 23days), (0.01, 20days)]

@info "=== Baroclinic adjustment: Ri parameter study ==="
@info "Fixed: f = $f, N² = $N², h = $h m, H = $H m"
@info "Waves: Uˢ = $(round(Uˢ, sigdigits=3)) m/s, ∂z_uˢ(0) = $(round(∂z_uˢ_surface, sigdigits=3)) s⁻¹"
@info "Ld = $(round(Ld, digits=1)) m"
@info "Grid: $Nx × $Ny × $Nz"
for (Ri, stop_time) in runs
    Ro = 1 / sqrt(Ri)
    Ps = f / (Ri * ∂z_uˢ_surface)
    @info @sprintf("  Ri = %g → Ro = %.3g, Ps = %.3g, stop_time = %s",
                   Ri, Ro, Ps, prettytime(stop_time))
end

for (Ri, stop_time) in runs
    run_simulation(case = :with_stokes; Ri, stop_time)
    run_simulation(case = :no_stokes; Ri, stop_time)
end

@info "All baroclinic adjustment simulations complete!"
