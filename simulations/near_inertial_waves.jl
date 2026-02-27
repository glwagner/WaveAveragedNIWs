# Near-inertial wave simulation with and without wave-averaged equations
#
# This script sets up two paired LES simulations following Wagner et al. (2021):
#
#   1. "growing_waves": Stokes drift grows in time, driving near-inertial waves
#      via the Craik-Leibovich vortex force and ∂t_uˢ. No surface stress.
#
#   2. "effective_stress": No Stokes drift. Instead, an equivalent time-dependent
#      surface stress τ(t) is applied that imparts the same total momentum as
#      the growing wave field.
#
# Reference:
#   Wagner, G. L., Chini, G. P., Ramadhan, A., Gallet, B., & Ferrari, R. (2021).
#   Near-inertial waves and turbulence driven by the growth of swell.
#   Journal of Physical Oceanography, 51(5), 1337-1351.

using Oceananigans
using Oceananigans.Units: minute, minutes, hours, hour
using Printf

# ============================
# Physical parameters
# ============================

const f = 1e-4       # Coriolis parameter [s⁻¹]
const N² = 1e-6      # Initial buoyancy gradient [s⁻²]
const g = Oceananigans.defaults.gravitational_acceleration

# Wave parameters (monochromatic deep-water surface gravity waves)
const wave_amplitude = 2.0    # m
const wave_wavelength = 100.0 # m
const wavenumber = 2π / wave_wavelength           # m⁻¹
const frequency = sqrt(g * wavenumber)             # s⁻¹
const vertical_scale = wave_wavelength / (4π)      # 1 / (2k) decay scale
const Uˢ_surface = wave_amplitude^2 * wavenumber * frequency  # surface Stokes drift [m s⁻¹]

# Growth time scale for swell
const T_growth = 4hours  # s

# Weak surface buoyancy flux to help spin up turbulence
const Qᵇ = 5e-10  # m² s⁻³ (destabilizing, positive = cooling in Oceananigans convention)

# Domain
const Lx = 128.0 # m
const Ly = 128.0 # m
const Lz = 64.0  # m

# ============================
# Resolution (low for CPU testing; increase for production)
# ============================
const Nx = 32
const Ny = 32
const Nz = 32

# ============================
# Stokes drift functions for the growing wave case
# ============================
#
# Stokes drift profile:  uˢ(z, t) = Uˢ_surface * exp(2k * z) * ramp(t)
# where ramp(t) = 1 - exp(-t² / (2T²))

@inline ramp(t)    = 1 - exp(-t^2 / (2 * T_growth^2))
@inline ∂t_ramp(t) = t / T_growth^2 * exp(-t^2 / (2 * T_growth^2))

# Vertical derivative of Stokes drift (needed for Craik-Leibovich vortex force)
@inline ∂z_uˢ_growing(z, t) = 2 * wavenumber * Uˢ_surface * exp(2 * wavenumber * z) * ramp(t)

# Time derivative of Stokes drift (needed for Lagrangian-mean momentum equation)
@inline ∂t_uˢ_growing(z, t) = Uˢ_surface * exp(2 * wavenumber * z) * ∂t_ramp(t)

# ============================
# Effective surface stress for the no-wave case
# ============================
#
# The vertically-integrated Stokes drift tendency equals the effective stress:
#   τ_eff(t) = ∫ ∂t_uˢ dz = Uˢ_surface / (2k) * ∂t_ramp(t)
#
# In Oceananigans, a negative surface flux drives a positive velocity.

@inline effective_stress(x, y, t) = -Uˢ_surface / (2 * wavenumber) * ∂t_ramp(t)

# ============================
# Sponge layer at the bottom
# ============================
# Relaxes velocities and buoyancy to their background state near z = -Lz

const sponge_width = 4.0   # m
const sponge_timescale = 60.0 # s

@inline μ(z) = exp(-(z + Lz)^2 / (2 * sponge_width^2)) / sponge_timescale

@inline sponge_u(x, y, z, t, u) = -μ(z) * u
@inline sponge_v(x, y, z, t, v) = -μ(z) * v
@inline sponge_w(x, y, z, t, w) = -μ(z) * w
@inline sponge_b(x, y, z, t, b) = -μ(z) * (b - N² * z)

# ============================
# Initial conditions
# ============================
# - Linear stratification: b(z) = N² z
# - Noise to seed turbulence, concentrated near the surface

const noise_amplitude = 1e-2  # Velocity noise amplitude [m s⁻¹]

uᵢ(x, y, z) = noise_amplitude * (2 * rand() - 1) * exp(z / 4)
vᵢ(x, y, z) = noise_amplitude * (2 * rand() - 1) * exp(z / 4)
wᵢ(x, y, z) = noise_amplitude * (2 * rand() - 1) * exp(z / 4)
bᵢ(x, y, z) = N² * z + 1e-4 * N² * Lz * (2 * rand() - 1) * exp(z / 4)

# ============================
# Build and run a simulation
# ============================

function run_simulation(; case, stop_time = 2 * 2π / f)

    grid = RectilinearGrid(CPU(),
                           size = (Nx, Ny, Nz),
                           extent = (Lx, Ly, Lz),
                           topology = (Periodic, Periodic, Bounded))

    coriolis = FPlane(; f)

    # Forcing: sponge layer
    u_forcing = Forcing(sponge_u; field_dependencies = :u)
    v_forcing = Forcing(sponge_v; field_dependencies = :v)
    w_forcing = Forcing(sponge_w; field_dependencies = :w)
    b_forcing = Forcing(sponge_b; field_dependencies = :b)

    forcing = (u = u_forcing, v = v_forcing, w = w_forcing, b = b_forcing)

    # Case-dependent setup
    if case == :growing_waves
        stokes_drift = UniformStokesDrift(∂z_uˢ = ∂z_uˢ_growing,
                                          ∂t_uˢ = ∂t_uˢ_growing)
        u_bcs = FieldBoundaryConditions()
    elseif case == :effective_stress
        stokes_drift = nothing
        u_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(effective_stress))
    else
        error("Unknown case: $case. Use :growing_waves or :effective_stress.")
    end

    b_bcs = FieldBoundaryConditions(top = FluxBoundaryCondition(Qᵇ),
                                    bottom = GradientBoundaryCondition(N²))

    boundary_conditions = (u = u_bcs, b = b_bcs)

    model = NonhydrostaticModel(grid; coriolis, forcing,
                                advection = WENO(),
                                tracers = :b,
                                buoyancy = BuoyancyTracer(),
                                closure = AnisotropicMinimumDissipation(),
                                stokes_drift,
                                boundary_conditions)

    # Set initial conditions
    set!(model, u = uᵢ, v = vᵢ, w = wᵢ, b = bᵢ)

    # Create simulation
    simulation = Simulation(model; Δt = 10.0, stop_time)
    conjure_time_step_wizard!(simulation, cfl = 0.5, max_Δt = 1minute)

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
                   schedule = TimeInterval(15minutes),
                   filename = "$(filename)_fields.jld2",
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

# ============================
# Run both cases
# ============================
#
# Run for 2 inertial periods: T_inertial = 2π/f

T_inertial = 2π / f
stop_time = 2 * T_inertial

@info "Inertial period: $(prettytime(T_inertial))"
@info "Stop time: $(prettytime(stop_time))"
@info "Grid: $Nx × $Ny × $Nz"
@info "Domain: $Lx × $Ly × $Lz m"

sim_waves = run_simulation(case = :growing_waves, stop_time = stop_time)
sim_stress = run_simulation(case = :effective_stress, stop_time = stop_time)

@info "Both simulations complete!"
