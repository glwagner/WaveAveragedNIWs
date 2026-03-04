# Langmuir turbulence LES after McWilliams, Sullivan & Moeng (1997)
#
# Two cases:
#   1. "langmuir": Wind stress + surface cooling + Stokes drift (CL vortex force)
#   2. "shear_convection": Wind stress + surface cooling, no Stokes drift
#
# Reference:
#   McWilliams, J. C., Sullivan, P. P., & Moeng, C.-H. (1997).
#   Langmuir turbulence in the ocean.
#   Journal of Fluid Mechanics, 334, 1-30.

using Oceananigans
using Oceananigans.Units: minute, minutes, hours
using Printf
using CUDA

# Physical parameters
f = 1e-4  # Coriolis parameter [s⁻¹]
g = Oceananigans.defaults.gravitational_acceleration

# Surface forcing (MSM97 values)
u★ = 6.1e-3               # friction velocity [m s⁻¹]
τx = -u★^2                # kinematic wind stress [m² s⁻²] (wind in +x)
Qᵇ = 2.307e-8             # destabilizing surface buoyancy flux [m² s⁻³]

# Stratification
N² = 1.936e-5              # buoyancy frequency squared [s⁻²]
initial_mixed_layer_depth = 33.0  # m

# Wave parameters (monochromatic deep-water surface gravity waves)
amplitude = 0.8            # m
wavelength = 60.0          # m
wavenumber = 2π / wavelength
frequency = sqrt(g * wavenumber)
vertical_scale = wavelength / (4π)  # 1/(2k) decay scale
Uˢ_surface = amplitude^2 * wavenumber * frequency

# Turbulent Langmuir number
La_t = sqrt(u★ / Uˢ_surface)

# Stokes drift
stokes_drift_parameters = (; Uˢ = Uˢ_surface, h = vertical_scale)
@inline ∂z_uˢ(z, t, p) = p.Uˢ / p.h * exp(z / p.h)

# Domain and grid
Lx = Ly = 256.0  # m
Lz = 128.0       # m
Nx = Ny = 128
Nz = 64

# Initial conditions
noise_amplitude = 1e-4  # m s⁻¹ (small perturbation, turbulence develops from forcing)

stratification(z) = z < -initial_mixed_layer_depth ? N² * z : N² * (-initial_mixed_layer_depth)

uᵢ(x, y, z) = noise_amplitude * (2 * rand() - 1) * exp(z / 4)
vᵢ(x, y, z) = noise_amplitude * (2 * rand() - 1) * exp(z / 4)
wᵢ(x, y, z) = noise_amplitude * (2 * rand() - 1) * exp(z / 4)
bᵢ(x, y, z) = stratification(z) + 1e-6 * (2 * rand() - 1) * exp(z / 4)

# Sponge layer at the bottom
sponge_width = 4.0    # m
sponge_timescale = 60.0 # s
sponge_params = (; Lz, δ = sponge_width, τ = sponge_timescale, N² = N²)

@inline μ(z, p) = exp(-(z + p.Lz)^2 / (2 * p.δ^2)) / p.τ

@inline sponge_u(x, y, z, t, u, p) = -μ(z, p) * u
@inline sponge_v(x, y, z, t, v, p) = -μ(z, p) * v
@inline sponge_w(x, y, z, t, w, p) = -μ(z, p) * w
@inline sponge_b(x, y, z, t, b, p) = -μ(z, p) * (b - p.N² * z)

function run_simulation(; case, stop_time = 12hours)

    grid = RectilinearGrid(GPU(),
                           size = (Nx, Ny, Nz),
                           extent = (Lx, Ly, Lz),
                           topology = (Periodic, Periodic, Bounded))

    coriolis = FPlane(; f)

    # Sponge layer forcing
    u_forcing = Forcing(sponge_u; field_dependencies = :u, parameters = sponge_params)
    v_forcing = Forcing(sponge_v; field_dependencies = :v, parameters = sponge_params)
    w_forcing = Forcing(sponge_w; field_dependencies = :w, parameters = sponge_params)
    b_forcing = Forcing(sponge_b; field_dependencies = :b, parameters = sponge_params)
    forcing = (u = u_forcing, v = v_forcing, w = w_forcing, b = b_forcing)

    # Case-dependent Stokes drift
    if case == :langmuir
        stokes_drift = UniformStokesDrift(∂z_uˢ = ∂z_uˢ, parameters = stokes_drift_parameters)
    elseif case == :shear_convection
        stokes_drift = nothing
    else
        error("Unknown case: $case. Use :langmuir or :shear_convection.")
    end

    # Boundary conditions: wind stress + surface cooling
    u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(τx))
    b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵇ),
                                    bottom = GradientBoundaryCondition(N²))

    boundary_conditions = (u = u_bcs, b = b_bcs)

    model = NonhydrostaticModel(grid; coriolis, forcing,
                                advection = WENO(order=9),
                                tracers = :b,
                                buoyancy = BuoyancyTracer(),
                                stokes_drift,
                                boundary_conditions)

    set!(model, u = uᵢ, v = vᵢ, w = wᵢ, b = bᵢ)

    simulation = Simulation(model; Δt = 1.0, stop_time)
    conjure_time_step_wizard!(simulation, cfl = 0.5, max_Δt = 30.0)

    # Progress callback
    function progress(sim)
        u, v, w = sim.model.velocities
        msg = @sprintf("Case: %s | i: %04d, t: %s, Δt: %s, max|u|: (%.1e, %.1e, %.1e) m/s, wall: %s\n",
                       case,
                       iteration(sim),
                       prettytime(time(sim)),
                       prettytime(sim.Δt),
                       maximum(abs, u), maximum(abs, v), maximum(abs, w),
                       prettytime(sim.run_wall_time))
        @info msg
    end

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(100))

    # Output
    filename = "langmuir_$(case)"
    u, v, w = model.velocities
    b = model.tracers.b

    # 2D slices for animation
    jslice = Ny ÷ 2
    w_xz = Field(w; indices = (:, jslice, :))
    u_xz = Field(u; indices = (:, jslice, :))
    b_xz = Field(b; indices = (:, jslice, :))
    w_xy = Field(w; indices = (:, :, Nz))
    u_xy = Field(u; indices = (:, :, Nz))

    simulation.output_writers[:slices] =
        JLD2Writer(model, (; w_xz, u_xz, b_xz, w_xy, u_xy),
                   schedule = TimeInterval(5minutes),
                   filename = "$(filename)_slices.jld2",
                   overwrite_existing = true)

    # Horizontally-averaged profiles
    U = Average(u, dims = (1, 2))
    V = Average(v, dims = (1, 2))
    B = Average(b, dims = (1, 2))
    wu = Average(w * u, dims = (1, 2))
    wv = Average(w * v, dims = (1, 2))
    wb = Average(w * b, dims = (1, 2))
    e = Average((u * u + v * v + w * w) / 2, dims = (1, 2))

    simulation.output_writers[:averages] =
        JLD2Writer(model, (; U, V, B, wu, wv, wb, e),
                   schedule = AveragedTimeInterval(5minutes, window = 2minutes),
                   filename = "$(filename)_averages.jld2",
                   overwrite_existing = true)

    @info "Running case: $case for $(prettytime(stop_time))..."
    @info "  La_t = $(round(La_t, digits=3))"
    @info "  u★ = $u★ m/s"
    @info "  Uˢ = $Uˢ_surface m/s"
    @info "  τx = $τx m²/s²"
    @info "  Qᵇ = $Qᵇ m²/s³"
    run!(simulation)

    @info "Case $case complete."
    return simulation
end

stop_time = 12hours

@info "=== MSM97 Langmuir turbulence ==="
@info "La_t = $(round(La_t, digits=3))"
@info "Grid: $Nx × $Ny × $Nz"
@info "Domain: $Lx × $Ly × $Lz m"

run_simulation(case = :langmuir, stop_time = stop_time)
run_simulation(case = :shear_convection, stop_time = stop_time)

@info "Both simulations complete!"
