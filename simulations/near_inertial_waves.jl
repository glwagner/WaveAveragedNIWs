# Near-inertial wave simulation with and without wave-averaged equations
#
# Three paired LES simulations based on the IVP in Section 4 of
# Wagner et al. (2021):
#
#   1. "wave_averaged": Steady Stokes drift with Craik-Leibovich equations.
#      Initial Lagrangian-mean velocity u_L = 1.1 uˢ(z).
#      Eulerian velocity is u_E = u_L - uˢ = 0.1 uˢ(z).
#
#   2. "no_waves_lagrangian": No Stokes drift.
#      Initial velocity u = 1.1 uˢ(z), matching the Lagrangian-mean
#      velocity of the wave-averaged case.
#
#   3. "no_waves_eulerian": No Stokes drift.
#      Initial velocity u = 0.1 uˢ(z), matching the Eulerian
#      velocity of the wave-averaged case.
#
# Reference:
#   Wagner, G. L., Chini, G. P., Ramadhan, A., Gallet, B., & Ferrari, R. (2021).
#   Near-inertial waves and turbulence driven by the growth of swell.
#   Journal of Physical Oceanography, 51(5), 1337-1351.

using Oceananigans
using Oceananigans.Units: minute, minutes, hours
using Printf
using CUDA

# Physical parameters
f = 1e-4       # Coriolis parameter [s⁻¹]
N² = 2e-5  # Initial buoyancy gradient [s⁻²]
g = Oceananigans.defaults.gravitational_acceleration

# Wave parameters (monochromatic deep-water surface gravity waves)
amplitude = 0.8  # m
wavelength = 60  # m
wavenumber = 2π / wavelength            # m⁻¹
frequency = sqrt(g * wavenumber)        # s⁻¹
vertical_scale = wavelength / (4π)      # 1 / (2k) decay scale
Uˢ_surface = amplitude^2 * wavenumber * frequency  # surface Stokes drift [m s⁻¹]

# Weak destabilizing surface buoyancy flux
Qᵇ = 1e-8  # m² s⁻³

# Surface kinematic momentum flux
τx = 0  # m² s⁻²

# Domain
Lx = Ly = 512 # m
Lz = 256 # m

# Resolution: 4m horizontal, 2m vertical
Nx = 128
Ny = 128
Nz = 128

# Stokes drift profile: uˢ(z) = Uˢ_surface * exp(z / vertical_scale)
uˢ(z) = Uˢ_surface * exp(z / vertical_scale)

# Vertical derivative of steady Stokes drift (for Craik-Leibovich vortex force)
stokes_drift_parameters = (; Uˢ = Uˢ_surface, h = vertical_scale)
@inline ∂z_uˢ(z, t, p) = p.Uˢ / p.h * exp(z / p.h)

# Sponge layer at the bottom: relaxes velocities and buoyancy near z = -Lz
sponge_width = 4.0    # m
sponge_timescale = 60.0 # s
sponge_params = (; Lz, δ = sponge_width, τ = sponge_timescale, N² = N²)

@inline μ(z, p) = exp(-(z + p.Lz)^2 / (2 * p.δ^2)) / p.τ

@inline sponge_u(x, y, z, t, u, p) = -μ(z, p) * u
@inline sponge_v(x, y, z, t, v, p) = -μ(z, p) * v
@inline sponge_w(x, y, z, t, w, p) = -μ(z, p) * w
@inline sponge_b(x, y, z, t, b, p) = -μ(z, p) * (b - p.N² * z)

# Initial conditions
noise_amplitude = 1e-2  # m s⁻¹
initial_mixed_layer_depth = 50

stratification(z) = z < -initial_mixed_layer_depth ? N² * z : N² * (-initial_mixed_layer_depth)

noise(z) = noise_amplitude * (2 * rand() - 1) * exp(z / 4)
vᵢ(x, y, z) = noise(z)
wᵢ(x, y, z) = noise(z)
bᵢ(x, y, z) = stratification(z) + 1e-1 * (2 * rand() - 1) * N² * Lz * exp(z / 4)

function run_simulation(; case, stop_time = 4 * 2π / f)

    grid = RectilinearGrid(GPU(),
                           size = (Nx, Ny, Nz),
                           extent = (Lx, Ly, Lz),
                           topology = (Periodic, Periodic, Bounded))

    coriolis = FPlane(; f)

    # Forcing: sponge layer
    u_forcing = Forcing(sponge_u; field_dependencies = :u, parameters = sponge_params)
    v_forcing = Forcing(sponge_v; field_dependencies = :v, parameters = sponge_params)
    w_forcing = Forcing(sponge_w; field_dependencies = :w, parameters = sponge_params)
    b_forcing = Forcing(sponge_b; field_dependencies = :b, parameters = sponge_params)

    forcing = (u = u_forcing, v = v_forcing, w = w_forcing, b = b_forcing)

    # Case-dependent setup
    if case == :wave_averaged
        stokes_drift = UniformStokesDrift(∂z_uˢ = ∂z_uˢ, parameters = stokes_drift_parameters)
        # Lagrangian-mean velocity: u_L = 1.1 uˢ (so Eulerian u_E = 0.1 uˢ)
        u_multiplier = 1.1
    elseif case == :no_waves_lagrangian
        stokes_drift = nothing
        # Same velocity as the Lagrangian-mean in wave_averaged case
        u_multiplier = 1.1
    elseif case == :no_waves_eulerian
        stokes_drift = nothing
        # Same velocity as the Eulerian velocity in wave_averaged case
        u_multiplier = 0.1
    else
        error("Unknown case: $case. Use :wave_averaged, :no_waves_lagrangian, or :no_waves_eulerian.")
    end

    uᵢ(x, y, z) = u_multiplier * uˢ(z) + noise(z)

    u_bcs = FieldBoundaryConditions()
    b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵇ),
                                    bottom = GradientBoundaryCondition(N²))

    boundary_conditions = (u = u_bcs, b = b_bcs)

    model = NonhydrostaticModel(grid; coriolis, forcing,
                                advection = WENO(order=9),
                                tracers = :b,
                                buoyancy = BuoyancyTracer(),
                                stokes_drift,
                                boundary_conditions)

    # Set initial conditions
    set!(model, u = uᵢ, v = vᵢ, w = wᵢ, b = bᵢ)

    # Create simulation
    simulation = Simulation(model; Δt = 10, stop_time)
    conjure_time_step_wizard!(simulation, cfl=0.7)

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

    simulation.callbacks[:progress] = Callback(progress, IterationInterval(20))

    # Output
    filename = "near_inertial_$(case)"

    u, v, w = model.velocities
    b = model.tracers.b

    # 3D snapshots (infrequent)
    simulation.output_writers[:fields] =
        JLD2Writer(model, (; u, v, w, b),
                   schedule = TimeInterval(2hours),
                   filename = "$(filename)_fields.jld2",
                   overwrite_existing = true)

    # 2D slices (frequent, for animation)
    w_xz = Field(w; indices = (:, 1, :))
    u_xz = Field(u; indices = (:, 1, :))
    b_xz = Field(b; indices = (:, 1, :))
    w_xy = Field(w; indices = (:, :, Nz))
    u_xy = Field(u; indices = (:, :, Nz))

    simulation.output_writers[:slices] =
        JLD2Writer(model, (; w_xz, u_xz, b_xz, w_xy, u_xy),
                   schedule = TimeInterval(5minutes),
                   filename = "$(filename)_slices.jld2",
                   overwrite_existing = true)

    # Horizontally-averaged profiles (more frequent)
    U = Average(u, dims = (1, 2))
    V = Average(v, dims = (1, 2))
    B = Average(b, dims = (1, 2))
    wu = Average(w * u, dims = (1, 2))
    wv = Average(w * v, dims = (1, 2))

    simulation.output_writers[:averages] =
        JLD2Writer(model, (; U, V, B, wu, wv),
                   schedule = AveragedTimeInterval(5minutes, window = 2minutes),
                   filename = "$(filename)_averages.jld2",
                   overwrite_existing = true)

    # Run!
    @info "Running case: $case for $(prettytime(stop_time))..."
    run!(simulation)

    @info "Case $case complete."
    return simulation
end

# Run for 2 inertial periods
T_inertial = 2π / f
stop_time = 2 * T_inertial

@info "Inertial period: $(prettytime(T_inertial))"
@info "Stop time: $(prettytime(stop_time))"
@info "Grid: $Nx × $Ny × $Nz"
@info "Domain: $Lx × $Ly × $Lz m"
@info "Surface Stokes drift: $Uˢ_surface m/s"

run_simulation(case = :wave_averaged, stop_time = stop_time)
run_simulation(case = :no_waves_lagrangian, stop_time = stop_time)
run_simulation(case = :no_waves_eulerian, stop_time = stop_time)

@info "All three simulations complete!"
