"""
AstroCalc Delta — Physics Engine
Sources:
  - NIST CODATA 2018 (constants)
  - Bate, Mueller & White, "Fundamentals of Astrodynamics" (orbital)
  - Sutton & Biblarz, "Rocket Propulsion Elements" 8th ed. (propulsion)
  - Murray & Dermott, "Solar System Dynamics" (planetary)
  - Misner, Thorne & Wheeler, "Gravitation" (relativity)
  - Carroll, "Spacetime and Geometry" (relativity/BH)
  - Griffiths, "Introduction to Electrodynamics" 4th ed. (EM)
  - Griffiths, "Introduction to Quantum Mechanics" 2nd ed. (QM)
  - Kittel, "Introduction to Solid State Physics" 8th ed. (solid state)
"""

import math

# ─────────────────────────────────────────────────────────────────────────────
# PHYSICAL CONSTANTS  (NIST CODATA 2018)
# ─────────────────────────────────────────────────────────────────────────────

class Constants:
    G         = 6.67430e-11       # Gravitational constant (m³ kg⁻¹ s⁻²)
    C         = 299_792_458.0     # Speed of light (m/s)
    H         = 6.62607015e-34    # Planck constant (J·s)
    H_BAR     = 1.054571817e-34   # Reduced Planck constant (J·s)
    K_B       = 1.380649e-23      # Boltzmann constant (J/K)
    E         = 1.602176634e-19   # Elementary charge (C)
    M_E       = 9.1093837015e-31  # Electron mass (kg)
    M_P       = 1.67262192369e-27 # Proton mass (kg)
    EPSILON_0 = 8.8541878128e-12  # Vacuum permittivity (F/m)
    MU_0      = 1.25663706212e-6  # Vacuum permeability (H/m)
    SIGMA     = 5.670374419e-8    # Stefan-Boltzmann (W m⁻² K⁻⁴)
    AU        = 1.495978707e11    # Astronomical unit (m)
    LY        = 9.4607304725808e15# Light year (m)
    PC        = 3.085677581e16    # Parsec (m)
    G_0       = 9.80665           # Standard gravity (m/s²)
    A_0       = 5.29177210903e-11 # Bohr radius (m)
    R_INF     = 1.0973731568539e7 # Rydberg constant (m⁻¹)
    ALPHA     = 7.2973525693e-3   # Fine structure constant

    # Earth
    M_EARTH   = 5.972168e24      # kg
    R_EARTH   = 6.371e6          # m
    MU_EARTH  = 3.986004418e14   # m³/s²
    RHO_ATM   = 1.225            # kg/m³ sea-level air density

    # Moon
    M_MOON    = 7.342e22         # kg
    R_MOON    = 1.7374e6         # m
    MU_MOON   = 4.9048695e12     # m³/s²
    A_MOON    = 3.844e8          # m (semi-major axis)

    # Sun
    M_SUN     = 1.989e30         # kg
    R_SUN     = 6.957e8          # m
    L_SUN     = 3.828e26         # W
    T_SUN     = 5778.0           # K (effective surface temp)

    # Mars
    M_MARS    = 6.4171e23        # kg
    R_MARS    = 3.3895e6         # m
    MU_MARS   = 4.282837e13      # m³/s²


# ─────────────────────────────────────────────────────────────────────────────
# DOMAIN 1: ORBITAL MECHANICS
# Ref: Bate, Mueller & White — Fundamentals of Astrodynamics (BMW 1971)
# ─────────────────────────────────────────────────────────────────────────────

class OrbitalMechanics:

    @staticmethod
    def orbital_velocity(altitude_m, mu=Constants.MU_EARTH, radius_body=Constants.R_EARTH):
        """Circular orbital velocity. v = sqrt(mu/r)  [BMW §2.3]"""
        return math.sqrt(mu / (radius_body + altitude_m))

    @staticmethod
    def orbital_period(altitude_m, mu=Constants.MU_EARTH, radius_body=Constants.R_EARTH):
        """Orbital period. T = 2pi*sqrt(r^3/mu)  [BMW §2.4]"""
        r = radius_body + altitude_m
        return 2 * math.pi * math.sqrt(r**3 / mu)

    @staticmethod
    def escape_velocity(altitude_m, mu=Constants.MU_EARTH, radius_body=Constants.R_EARTH):
        """Escape velocity. v_esc = sqrt(2mu/r)  [BMW §2.5]"""
        return math.sqrt(2 * mu / (radius_body + altitude_m))

    @staticmethod
    def hohmann_delta_v(r1, r2, mu=Constants.MU_EARTH):
        """Hohmann transfer delta-v. r1,r2 in metres.  [BMW §6.3]"""
        v1  = math.sqrt(mu / r1)
        v2  = math.sqrt(mu / r2)
        vt1 = math.sqrt(mu * (2/r1 - 2/(r1+r2)))
        vt2 = math.sqrt(mu * (2/r2 - 2/(r1+r2)))
        dv1 = abs(vt1 - v1)
        dv2 = abs(v2 - vt2)
        return {"dv1_ms": round(dv1,4), "dv2_ms": round(dv2,4),
                "total_dv_ms": round(dv1+dv2,4)}

    @staticmethod
    def vis_viva(r, a, mu=Constants.MU_EARTH):
        """Vis-viva. v = sqrt(mu*(2/r - 1/a))  [BMW §2.8]"""
        return math.sqrt(mu * (2/r - 1/a))

    @staticmethod
    def eccentricity(r_apoapsis, r_periapsis):
        """Eccentricity. e = (ra-rp)/(ra+rp)  [BMW §1.5]"""
        return (r_apoapsis - r_periapsis) / (r_apoapsis + r_periapsis)

    @staticmethod
    def inclination_change_dv(v_orbit, delta_i_deg):
        """Pure inclination change. dv = 2v*sin(di/2)  [BMW §6.6]"""
        return 2 * v_orbit * math.sin(math.radians(delta_i_deg) / 2)

    @staticmethod
    def geostationary_altitude(mu=Constants.MU_EARTH, omega=7.2921150e-5):
        """GEO altitude. r = (mu/omega^2)^(1/3) - R_E  [BMW §2.4]"""
        return (mu / omega**2) ** (1/3) - Constants.R_EARTH

    @staticmethod
    def sphere_of_influence(a, m_body, m_central=Constants.M_SUN):
        """SOI radius. r_SOI = a*(m/M)^(2/5)  [BMW §7.3]"""
        return a * (m_body / m_central) ** (2/5)

    @staticmethod
    def orbital_energy(r, v, mu=Constants.MU_EARTH):
        """Specific orbital energy. eps = v^2/2 - mu/r  [BMW §2.6]"""
        return v**2 / 2 - mu / r

    @staticmethod
    def semi_major_axis_from_period(T, mu=Constants.MU_EARTH):
        """SMA from period. a = (mu*(T/2pi)^2)^(1/3)  [BMW §2.4]"""
        return (mu * (T / (2*math.pi))**2) ** (1/3)

    # ── NEW ──────────────────────────────────────────────────────────────────

    @staticmethod
    def bi_elliptic_delta_v(r1, r2, r_transfer, mu=Constants.MU_EARTH):
        """
        Bi-elliptic transfer total delta-v. More efficient than Hohmann when r2/r1 > 11.94.
        r_transfer: apoapsis of intermediate ellipse (m).  [BMW §6.4]
        """
        v1       = math.sqrt(mu / r1)
        vt1_peri = math.sqrt(mu * (2/r1 - 2/(r1 + r_transfer)))
        vt1_apo  = math.sqrt(mu * (2/r_transfer - 2/(r1 + r_transfer)))
        vt2_apo  = math.sqrt(mu * (2/r_transfer - 2/(r_transfer + r2)))
        vt2_peri = math.sqrt(mu * (2/r2 - 2/(r_transfer + r2)))
        v2       = math.sqrt(mu / r2)
        dv1  = abs(vt1_peri - v1)
        dv2  = abs(vt2_apo - vt1_apo)
        dv3  = abs(v2 - vt2_peri)
        total = dv1 + dv2 + dv3
        return {"dv1_ms": round(dv1,4), "dv2_ms": round(dv2,4),
                "dv3_ms": round(dv3,4), "total_dv_ms": round(total,4)}

    @staticmethod
    def orbital_decay_timescale(altitude_m, Cd=2.2, A_m2=10.0, mass_kg=1000.0):
        """
        Characteristic drag decay timescale using exponential atmosphere.
        Scale height H = 8500 m, rho_0 = 1.225 kg/m^3.  [Vallado §9.6]
        Returns timescale in seconds.
        """
        v = math.sqrt(Constants.MU_EARTH / (Constants.R_EARTH + altitude_m))
        rho = Constants.RHO_ATM * math.exp(-altitude_m / 8500.0)
        a_drag = 0.5 * rho * v**2 * (Cd * A_m2 / mass_kg)
        if a_drag == 0:
            return float('inf')
        return v / (2 * a_drag)

    @staticmethod
    def repeating_ground_track(k_orbits, m_days, mu=Constants.MU_EARTH):
        """
        SMA for a repeating ground track: satellite completes k orbits in m sidereal days.
        Sidereal day = 86164.1 s.  [Wertz §6.4]
        """
        T_sat = m_days * 86164.1 / k_orbits
        return (mu * (T_sat / (2 * math.pi))**2) ** (1/3)

    @staticmethod
    def delta_v_plane_change_combined(v1, v2, delta_i_deg):
        """
        Combined plane change + speed change delta-v.
        dv = sqrt(v1^2 + v2^2 - 2*v1*v2*cos(di))  [BMW §6.6]
        """
        return math.sqrt(v1**2 + v2**2 - 2*v1*v2*math.cos(math.radians(delta_i_deg)))

    @staticmethod
    def nodal_regression_rate(altitude_m, inclination_deg, e=0.0):
        """
        J2 RAAN drift rate. dOmega/dt in deg/day.
        J2 = 1.08263e-3.  [Wertz §6.4; Vallado §9.3]
        """
        J2 = 1.08263e-3
        r  = Constants.R_EARTH + altitude_m
        n  = math.sqrt(Constants.MU_EARTH / r**3)
        rate_rad_s = (-3/2) * n * J2 * (Constants.R_EARTH/r)**2 * math.cos(math.radians(inclination_deg)) / (1-e**2)**2
        return math.degrees(rate_rad_s) * 86400


# ─────────────────────────────────────────────────────────────────────────────
# DOMAIN 2: ROCKET PROPULSION
# Ref: Sutton & Biblarz — Rocket Propulsion Elements, 8th ed. (S&B)
# ─────────────────────────────────────────────────────────────────────────────

class RocketPropulsion:

    @staticmethod
    def tsiolkovsky(isp, m0, mf):
        """Tsiolkovsky. dv = Isp*g0*ln(m0/mf)  [S&B §2.4]"""
        if mf <= 0 or m0 <= mf:
            raise ValueError("Masses invalid: m0 must be > mf > 0")
        return isp * Constants.G_0 * math.log(m0 / mf)

    @staticmethod
    def mass_ratio(delta_v, isp):
        """Mass ratio. m0/mf = e^(dv/ve)  [S&B §2.4]"""
        return math.exp(delta_v / (isp * Constants.G_0))

    @staticmethod
    def propellant_fraction(delta_v, isp):
        """Propellant fraction. 1 - e^(-dv/ve)  [S&B §2.4]"""
        return 1 - math.exp(-delta_v / (isp * Constants.G_0))

    @staticmethod
    def thrust(mdot, isp):
        """Thrust. F = mdot*Isp*g0  [S&B §3.1]"""
        return mdot * isp * Constants.G_0

    @staticmethod
    def specific_impulse(thrust_N, mdot_kgs):
        """Isp. Isp = F/(mdot*g0)  [S&B §3.1]"""
        return thrust_N / (mdot_kgs * Constants.G_0)

    @staticmethod
    def exhaust_velocity(isp):
        """Exhaust velocity. ve = Isp*g0  [S&B §2.3]"""
        return isp * Constants.G_0

    @staticmethod
    def burn_time(m0, mf, mdot):
        """Burn time. t = (m0-mf)/mdot  [S&B §2.4]"""
        return (m0 - mf) / mdot

    @staticmethod
    def twr(thrust_N, mass_kg, g=Constants.G_0):
        """Thrust-to-weight ratio.  [S&B §17.2]"""
        return thrust_N / (mass_kg * g)

    @staticmethod
    def gravity_loss_estimate(twr, burn_time_s):
        """Gravity loss. dv_loss = g0*t*(1-1/TWR)  [S&B §4.3]"""
        return Constants.G_0 * burn_time_s * (1 - 1/twr) if twr > 1 else float('inf')

    @staticmethod
    def staged_delta_v(stages):
        """Multi-stage delta-v sum.  [S&B §4.2]"""
        return sum(RocketPropulsion.tsiolkovsky(s["isp"], s["m0"], s["mf"]) for s in stages)

    # ── NEW ──────────────────────────────────────────────────────────────────

    @staticmethod
    def landing_burn_delta_v(mass_dry_kg, mass_prop_landing_kg, isp):
        """
        Delta-v available from landing propellant load (e.g. Falcon 9 boostback/landing).
        Uses Tsiolkovsky with landing prop only.  [S&B §2.4]
        """
        m0 = mass_dry_kg + mass_prop_landing_kg
        mf = mass_dry_kg
        return RocketPropulsion.tsiolkovsky(isp, m0, mf)

    @staticmethod
    def nozzle_exit_velocity(T_chamber_K, gamma, M_mol_kgmol, p_e_Pa, p_c_Pa):
        """
        Ideal nozzle exit velocity.
        ve = sqrt(2*gamma/(gamma-1) * R*T_c/M * [1-(pe/pc)^((gamma-1)/gamma)])
        [S&B §3.3]
        """
        R_u = 8.314462  # J/(mol·K)
        g   = gamma
        arg = (2*g/(g-1)) * (R_u * T_chamber_K / M_mol_kgmol) * (1 - (p_e_Pa/p_c_Pa)**((g-1)/g))
        return math.sqrt(max(arg, 0))

    @staticmethod
    def reuse_propellant_margin(m0_kg, isp, dv_ascent, dv_reentry, dv_landing, mf_dry_kg):
        """
        Remaining propellant after ascent + reentry + landing burns.
        Positive margin = vehicle is reusable on that mission.  [S&B §2.4]
        """
        ve = isp * Constants.G_0
        m1 = m0_kg    * math.exp(-dv_ascent  / ve)
        m2 = m1       * math.exp(-dv_reentry / ve)
        m3 = m2       * math.exp(-dv_landing / ve)
        return m3 - mf_dry_kg

    @staticmethod
    def optimal_staging_payload_fraction(dv_total, isp, n_stages, eps=0.1):
        """
        Payload fraction for optimal equal-dv n-stage rocket.
        eps: structural mass fraction per stage.  [S&B §4.2]
        """
        dv_stage = dv_total / n_stages
        ve = isp * Constants.G_0
        mr = math.exp(dv_stage / ve)
        pf_stage = (1 - eps) / mr - eps
        return {"mass_ratio_per_stage": round(mr, 4),
                "payload_frac_per_stage": round(max(pf_stage, 0), 6),
                "total_payload_frac": round(max(pf_stage**n_stages, 0), 8)}

    @staticmethod
    def chamber_pressure_thrust_coefficient(gamma, p_e, p_c):
        """
        Ideal thrust coefficient C_F (dimensionless).
        C_F = sqrt(2g^2/(g-1)*(2/(g+1))^((g+1)/(g-1))*(1-(pe/pc)^((g-1)/g)))
        [S&B §3.4]
        """
        g = gamma
        term1 = 2*g**2 / (g-1)
        term2 = (2/(g+1)) ** ((g+1)/(g-1))
        term3 = 1 - (p_e/p_c)**((g-1)/g)
        return math.sqrt(term1 * term2 * term3)


# ─────────────────────────────────────────────────────────────────────────────
# DOMAIN 3: PLANETARY SCIENCE & KEPLER
# Ref: Murray & Dermott — Solar System Dynamics (M&D 1999)
# ─────────────────────────────────────────────────────────────────────────────

class PlanetaryScience:

    @staticmethod
    def surface_gravity(mass_kg, radius_m):
        """Surface gravity. g = GM/R^2  [M&D §1.2]"""
        return Constants.G * mass_kg / radius_m**2

    @staticmethod
    def kepler_third_law_period(a_m):
        """Heliocentric period. T = 2pi*sqrt(a^3/mu_sun)  [M&D §2.2]"""
        return 2 * math.pi * math.sqrt(a_m**3 / (Constants.G * Constants.M_SUN))

    @staticmethod
    def kepler_third_law_sma(T_s):
        """SMA from heliocentric period.  [M&D §2.2]"""
        return (Constants.G * Constants.M_SUN * (T_s/(2*math.pi))**2) ** (1/3)

    @staticmethod
    def roche_limit(R_primary, rho_primary, rho_satellite):
        """Rigid-body Roche limit. d = R_p*(2*rho_p/rho_s)^(1/3)  [M&D §4.4]"""
        return R_primary * (2 * rho_primary / rho_satellite) ** (1/3)

    @staticmethod
    def hill_sphere(a, e, m_planet, m_star=Constants.M_SUN):
        """Hill sphere. r_H = a*(1-e)*(m/3M)^(1/3)  [M&D §3.3]"""
        return a * (1-e) * (m_planet / (3*m_star)) ** (1/3)

    @staticmethod
    def equilibrium_temperature(albedo, distance_au):
        """Equilibrium temp. T_eq = T_sun*sqrt(R_sun/2d)*(1-A)^0.25  [M&D §9.2]"""
        d = distance_au * Constants.AU
        return Constants.T_SUN * math.sqrt(Constants.R_SUN / (2*d)) * (1-albedo)**0.25

    @staticmethod
    def synodic_period(T1, T2):
        """Synodic period. 1/T_syn = |1/T1 - 1/T2|  [M&D §2.3]"""
        return 1 / abs(1/T1 - 1/T2)

    @staticmethod
    def planet_density(mass_kg, radius_m):
        """Mean density. rho = M/(4/3*pi*R^3)  [M&D §1.1]"""
        return mass_kg / ((4/3) * math.pi * radius_m**3)

    @staticmethod
    def tidal_force(m_tidal, r_center, delta_r):
        """Tidal acceleration. a = 2GM*dr/r^3  [M&D §4.1]"""
        return 2 * Constants.G * m_tidal * delta_r / r_center**3

    @staticmethod
    def angular_diameter(diameter_m, distance_m):
        """Angular diameter in arcseconds.  [M&D §1.3]"""
        return math.degrees(2 * math.atan(diameter_m / (2*distance_m))) * 3600

    # ── NEW ──────────────────────────────────────────────────────────────────

    @staticmethod
    def lunar_transfer_delta_v(parking_altitude_m=200e3):
        """
        Trans-lunar injection delta-v from LEO parking orbit.
        Patched-conic approximation.  [BMW §7.4]
        """
        r_park = Constants.R_EARTH + parking_altitude_m
        a_tli  = (r_park + Constants.A_MOON) / 2
        v_park = math.sqrt(Constants.MU_EARTH / r_park)
        v_tli  = math.sqrt(Constants.MU_EARTH * (2/r_park - 1/a_tli))
        return v_tli - v_park

    @staticmethod
    def mars_surface_weight(mass_kg):
        """Weight on Mars. W = m*g_Mars  [M&D §1.2]"""
        g_mars = Constants.G * Constants.M_MARS / Constants.R_MARS**2
        return mass_kg * g_mars

    @staticmethod
    def asteroid_impact_energy(mass_kg, velocity_ms):
        """
        Kinetic energy of asteroid impact. E = 0.5*m*v^2
        [Collins et al., MAPS 2005; Holsapple 1993]
        """
        E_J  = 0.5 * mass_kg * velocity_ms**2
        E_Mt = E_J / 4.184e15
        return {"energy_J": E_J, "energy_Mt_TNT": E_Mt}

    @staticmethod
    def lunar_escape_velocity():
        """Escape velocity from lunar surface. v = sqrt(2*mu_Moon/R_Moon)  [BMW §2.5]"""
        return math.sqrt(2 * Constants.MU_MOON / Constants.R_MOON)

    @staticmethod
    def escape_velocity_general(mass_kg, radius_m):
        """General escape velocity. v = sqrt(2GM/R)  [BMW §2.5]"""
        return math.sqrt(2 * Constants.G * mass_kg / radius_m)

    @staticmethod
    def orbital_resonance_ratio(T1_days, T2_days):
        """Near-integer resonance ratio T2/T1.  [M&D §8.1]"""
        return T2_days / T1_days


# ─────────────────────────────────────────────────────────────────────────────
# DOMAIN 4: RELATIVITY & BLACK HOLES
# Ref: MTW — Gravitation (1973); Carroll — Spacetime & Geometry (2004)
# ─────────────────────────────────────────────────────────────────────────────

class Relativity:

    @staticmethod
    def lorentz_factor(v):
        """Lorentz factor. gamma = 1/sqrt(1-v^2/c^2)  [MTW §2.1]"""
        beta = v / Constants.C
        if beta >= 1:
            raise ValueError("v must be < c")
        return 1 / math.sqrt(1 - beta**2)

    @staticmethod
    def time_dilation(t0, v):
        """Time dilation. t = gamma*t0  [MTW §2.2]"""
        return t0 * Relativity.lorentz_factor(v)

    @staticmethod
    def length_contraction(L0, v):
        """Length contraction. L = L0/gamma  [MTW §2.2]"""
        return L0 / Relativity.lorentz_factor(v)

    @staticmethod
    def relativistic_momentum(m, v):
        """Relativistic momentum. p = gamma*m*v  [MTW §2.3]"""
        return Relativity.lorentz_factor(v) * m * v

    @staticmethod
    def relativistic_energy(m, v):
        """Total energy. E = gamma*m*c^2  [MTW §2.3]"""
        return Relativity.lorentz_factor(v) * m * Constants.C**2

    @staticmethod
    def rest_energy(m):
        """Rest energy. E0 = m*c^2  [MTW §2.3]"""
        return m * Constants.C**2

    @staticmethod
    def relativistic_kinetic_energy(m, v):
        """Relativistic KE. K = (gamma-1)*m*c^2  [MTW §2.3]"""
        return (Relativity.lorentz_factor(v) - 1) * m * Constants.C**2

    @staticmethod
    def velocity_addition(u, v):
        """Velocity addition. w = (u+v)/(1+uv/c^2)  [MTW §2.4]"""
        return (u + v) / (1 + u*v / Constants.C**2)

    @staticmethod
    def schwarzschild_radius(M):
        """Schwarzschild radius. rs = 2GM/c^2  [MTW §23.1]"""
        return 2 * Constants.G * M / Constants.C**2

    @staticmethod
    def photon_sphere_radius(M):
        """Photon sphere. r_ph = 3GM/c^2  [MTW §25.5]"""
        return 3 * Constants.G * M / Constants.C**2

    @staticmethod
    def innermost_stable_orbit(M):
        """ISCO (Schwarzschild). r_ISCO = 6GM/c^2  [MTW §25.5]"""
        return 6 * Constants.G * M / Constants.C**2

    @staticmethod
    def kerr_event_horizon(M, a):
        """Kerr outer horizon. r+ = GM/c^2 + sqrt((GM/c^2)^2 - a^2)  [Carroll §9.6]"""
        rg   = Constants.G * M / Constants.C**2
        disc = rg**2 - a**2
        return None if disc < 0 else rg + math.sqrt(disc)

    @staticmethod
    def hawking_temperature(M):
        """Hawking temperature. T_H = hbar*c^3/(8pi*G*M*k_B)  [Hawking 1975]"""
        return Constants.H_BAR * Constants.C**3 / (8 * math.pi * Constants.G * M * Constants.K_B)

    @staticmethod
    def hawking_luminosity(M):
        """Hawking power. P = hbar*c^6/(15360*pi*G^2*M^2)  [Page 1976]"""
        return Constants.H_BAR * Constants.C**6 / (15360 * math.pi * Constants.G**2 * M**2)

    @staticmethod
    def bh_evaporation_time(M):
        """BH evaporation time. t = 5120*pi*G^2*M^3/(hbar*c^4)  [Page 1976]"""
        return 5120 * math.pi * Constants.G**2 * M**3 / (Constants.H_BAR * Constants.C**4)

    @staticmethod
    def gravitational_redshift(r_emit, M):
        """Gravitational redshift. z+1 = 1/sqrt(1-rs/r)  [MTW §25.2]"""
        rs = Relativity.schwarzschild_radius(M)
        if r_emit <= rs:
            return float('inf')
        return 1 / math.sqrt(1 - rs/r_emit)

    @staticmethod
    def light_deflection(M, b):
        """GR light deflection. alpha = 4GM/(c^2*b)  [MTW §40.3]"""
        return 4 * Constants.G * M / (Constants.C**2 * b)

    @staticmethod
    def gravitational_time_dilation(r, M):
        """Gravitational time dilation. sqrt(1-rs/r)  [MTW §25.2]"""
        rs = Relativity.schwarzschild_radius(M)
        if r <= rs:
            return 0.0
        return math.sqrt(1 - rs/r)

    # ── NEW ──────────────────────────────────────────────────────────────────

    @staticmethod
    def gravitational_wave_strain(M_chirp_kg, distance_m, frequency_hz):
        """
        GW strain amplitude h from compact binary inspiral (quadrupole formula).
        h ~ (4/r)*(G*M_chirp/c^2)^(5/3)*(pi*f/c)^(2/3)
        [Abbott et al. PRL 2016; Carroll §7.5]
        """
        G, C = Constants.G, Constants.C
        return (4/distance_m) * (G*M_chirp_kg/C**2)**(5/3) * (math.pi*frequency_hz/C)**(2/3)

    @staticmethod
    def lense_thirring_precession(r, J):
        """
        Lense-Thirring (frame dragging) precession rate.
        Omega_LT = 2GJ/(c^2*r^3)  rad/s
        J: angular momentum of central body (kg*m^2/s).  [MTW §40.7]
        """
        return 2 * Constants.G * J / (Constants.C**2 * r**3)

    @staticmethod
    def relativistic_precession_per_orbit(a_m, e, M_central_kg):
        """
        GR periapsis precession per orbit (Schwarzschild).
        Delta_phi = 6pi*G*M / (c^2*a*(1-e^2))  radians per orbit  [MTW §40.5]
        Returns arcseconds per orbit.
        """
        delta_phi = 6 * math.pi * Constants.G * M_central_kg / (Constants.C**2 * a_m * (1-e**2))
        return math.degrees(delta_phi) * 3600

    @staticmethod
    def photon_redshift_at_isco(M):
        """
        Gravitational redshift of photon emitted at ISCO.
        z = 1/sqrt(1 - rs/r_ISCO) - 1  [Carroll §9.5]
        For Schwarzschild: r_ISCO = 6GM/c^2, rs = 2GM/c^2, so 1-rs/r_ISCO = 2/3.
        """
        return 1 / math.sqrt(2/3) - 1

    @staticmethod
    def penrose_max_efficiency(a_star):
        """
        Max energy extraction efficiency from Kerr BH via Penrose process.
        Interpolated between Schwarzschild (1-1/sqrt(6)) and max Kerr (1-1/sqrt(3)).
        a_star: dimensionless spin [0,1].  [Carroll §9.7]
        """
        if not 0 <= a_star <= 1:
            raise ValueError("Spin a* must be in [0,1]")
        eta_s = 1 - 1/math.sqrt(6)
        eta_k = 1 - 1/math.sqrt(3)
        return eta_s + a_star * (eta_k - eta_s)


# ─────────────────────────────────────────────────────────────────────────────
# DOMAIN 5: ELECTROMAGNETISM
# Ref: Griffiths — Introduction to Electrodynamics, 4th ed.
# ─────────────────────────────────────────────────────────────────────────────

class ElectroMagnetism:

    @staticmethod
    def coulomb_force(q1, q2, r):
        """Coulomb's law. F = k*q1*q2/r^2  [Griffiths §2.1]"""
        return q1 * q2 / (4 * math.pi * Constants.EPSILON_0 * r**2)

    @staticmethod
    def electric_field_point(q, r):
        """E-field. E = kq/r^2  [Griffiths §2.1]"""
        return q / (4 * math.pi * Constants.EPSILON_0 * r**2)

    @staticmethod
    def electric_potential(q, r):
        """Electric potential. V = kq/r  [Griffiths §2.3]"""
        return q / (4 * math.pi * Constants.EPSILON_0 * r)

    @staticmethod
    def magnetic_force(q, v, B, theta_deg=90):
        """Lorentz force. F = qvB*sin(theta)  [Griffiths §5.1]"""
        return abs(q) * v * B * math.sin(math.radians(theta_deg))

    @staticmethod
    def cyclotron_radius(m, v, q, B):
        """Cyclotron radius. r = mv/(qB)  [Griffiths §5.2]"""
        return m * v / (abs(q) * B)

    @staticmethod
    def cyclotron_frequency(q, B, m):
        """Cyclotron frequency. f = qB/(2pi*m)  [Griffiths §5.2]"""
        return abs(q) * B / (2 * math.pi * m)

    @staticmethod
    def biot_savart_wire(I, r):
        """B-field from wire. B = mu0*I/(2pi*r)  [Griffiths §5.3]"""
        return Constants.MU_0 * I / (2 * math.pi * r)

    @staticmethod
    def poynting_vector(E, B):
        """Poynting vector. S = E*B/mu0  [Griffiths §8.1]"""
        return E * B / Constants.MU_0

    @staticmethod
    def skin_depth(rho, f, mu_r=1):
        """Skin depth. delta = sqrt(2*rho/(omega*mu))  [Griffiths §9.4]"""
        return math.sqrt(2 * rho / (2 * math.pi * f * mu_r * Constants.MU_0))

    @staticmethod
    def larmor_radiation_power(q, a):
        """Larmor radiation. P = q^2*a^2/(6pi*eps0*c^3)  [Griffiths §11.2]"""
        return q**2 * a**2 / (6 * math.pi * Constants.EPSILON_0 * Constants.C**3)

    @staticmethod
    def em_wave_speed(epsilon_r=1, mu_r=1):
        """EM wave speed. v = c/sqrt(eps_r*mu_r)  [Griffiths §9.3]"""
        return Constants.C / math.sqrt(epsilon_r * mu_r)

    # ── NEW ──────────────────────────────────────────────────────────────────

    @staticmethod
    def plasma_frequency(n_e):
        """
        Plasma angular frequency. omega_p = sqrt(n_e*e^2/(eps0*m_e))  rad/s
        n_e: electron number density (m^-3).
        [Chen, Plasma Physics §4.2; Griffiths §9.4]
        """
        return math.sqrt(n_e * Constants.E**2 / (Constants.EPSILON_0 * Constants.M_E))

    @staticmethod
    def magnetic_pressure(B):
        """
        Magnetic pressure. P_B = B^2/(2*mu0)  Pa
        [Griffiths §8.2; Freidberg, Plasma Physics §3.1]
        """
        return B**2 / (2 * Constants.MU_0)

    @staticmethod
    def synchrotron_radiation_power(gamma, B):
        """
        Synchrotron radiation power (relativistic electron).
        P = (q^4*B^2*gamma^2)/(6pi*eps0*m_e^2*c)
        [Jackson §14.4; Rybicki & Lightman §6.1]
        """
        q = Constants.E
        return (q**4 * B**2 * gamma**2) / (6 * math.pi * Constants.EPSILON_0 * Constants.M_E**2 * Constants.C)

    @staticmethod
    def alfven_wave_speed(B, rho_mass):
        """
        Alfven wave speed. v_A = B/sqrt(mu0*rho)  m/s
        [Chen §4.5; Priest, Solar MHD §2.3]
        """
        return B / math.sqrt(Constants.MU_0 * rho_mass)

    @staticmethod
    def hall_parameter(B, m, nu_collision):
        """
        Hall parameter (electron magnetization). beta_H = eB/(m*nu_c)
        [Chen §2.3; NRL Plasma Formulary]
        """
        return Constants.E * abs(B) / (m * nu_collision)


# ─────────────────────────────────────────────────────────────────────────────
# DOMAIN 6: QUANTUM MECHANICS & SOLID STATE
# Ref: Griffiths QM; Kittel Solid State; BCS 1957
# ─────────────────────────────────────────────────────────────────────────────

class QuantumMechanics:

    @staticmethod
    def photon_energy(freq):
        """E = h*f  [Griffiths §1.1]"""
        return Constants.H * freq

    @staticmethod
    def photon_energy_from_wavelength(wavelength_m):
        """E = hc/lambda  [Griffiths §1.1]"""
        return Constants.H * Constants.C / wavelength_m

    @staticmethod
    def photon_frequency(wavelength_m):
        """f = c/lambda"""
        return Constants.C / wavelength_m

    @staticmethod
    def de_broglie_wavelength(m, v):
        """lambda = h/(mv)  [Griffiths §1.2]"""
        return Constants.H / (m * v)

    @staticmethod
    def de_broglie_wavelength_relativistic(m, v):
        """Relativistic de Broglie. lambda = h/(gamma*m*v)  [Griffiths §1.2]"""
        return Constants.H / Relativity.relativistic_momentum(m, v)

    @staticmethod
    def photoelectric_max_ke(freq, work_function_eV):
        """KE = hf - phi  [Griffiths §1.1]"""
        return max(0.0, Constants.H * freq - work_function_eV * Constants.E)

    @staticmethod
    def hydrogen_energy_level(n):
        """E_n = -13.6/n^2 eV  [Griffiths §4.2]"""
        return -13.6 / n**2

    @staticmethod
    def hydrogen_transition_wavelength(n1, n2):
        """1/lambda = R_H*(1/n1^2 - 1/n2^2)  [Griffiths §4.2]"""
        inv = Constants.R_INF * (1/n1**2 - 1/n2**2)
        return None if inv <= 0 else 1 / inv

    @staticmethod
    def uncertainty_position_momentum(delta_p):
        """dx >= hbar/(2*dp)  [Griffiths §3.5]"""
        return Constants.H_BAR / (2 * delta_p)

    @staticmethod
    def uncertainty_energy_time(delta_t):
        """dE >= hbar/(2*dt)  [Griffiths §3.5]"""
        return Constants.H_BAR / (2 * delta_t)

    @staticmethod
    def compton_wavelength_shift(theta_deg):
        """Dlambda = (h/m_e*c)*(1-cos(theta))  [Griffiths §1.3]"""
        return (Constants.H / (Constants.M_E * Constants.C)) * (1 - math.cos(math.radians(theta_deg)))

    @staticmethod
    def stefan_boltzmann(T):
        """P = sigma*T^4  [Griffiths §1.1]"""
        return Constants.SIGMA * T**4

    @staticmethod
    def wien_peak_wavelength(T):
        """lambda_max = b/T  [Griffiths §1.1]"""
        return 2.897771955e-3 / T

    @staticmethod
    def bohr_radius_n(n):
        """r_n = n^2*a0  [Griffiths §4.2]"""
        return n**2 * Constants.A_0

    # ── NEW ──────────────────────────────────────────────────────────────────

    @staticmethod
    def fermi_energy(n_e):
        """
        Free electron Fermi energy.
        E_F = (hbar^2/2m_e)*(3pi^2*n)^(2/3)  [Kittel §6.1]
        n_e: electron density (m^-3). Returns Joules.
        """
        return (Constants.H_BAR**2 / (2 * Constants.M_E)) * (3 * math.pi**2 * n_e)**(2/3)

    @staticmethod
    def cooper_pair_binding_energy(T_c):
        """
        BCS superconducting gap (Cooper pair binding energy).
        2*Delta = 3.528*k_B*T_c  [BCS 1957; Kittel §12.5]
        Returns gap energy in Joules.
        """
        return 3.528 * Constants.K_B * T_c

    @staticmethod
    def josephson_frequency(V):
        """
        AC Josephson frequency. f = 2eV/h  [Josephson 1962; Kittel §12.6]
        V: voltage across junction (V). Returns Hz.
        """
        return 2 * Constants.E * V / Constants.H

    @staticmethod
    def quantum_tunneling_probability(m, V0_J, E_J, L_m):
        """
        WKB tunneling probability through rectangular barrier.
        T ~ exp(-2*kappa*L),  kappa = sqrt(2m(V0-E))/hbar  [Griffiths §8.1]
        Returns transmission coefficient in [0,1].
        """
        if E_J >= V0_J:
            return 1.0
        kappa = math.sqrt(2 * m * (V0_J - E_J)) / Constants.H_BAR
        return math.exp(-2 * kappa * L_m)

    @staticmethod
    def debye_temperature(v_sound, n_density):
        """
        Debye temperature. Theta_D = (hbar/k_B)*v_s*(6pi^2*n)^(1/3)
        [Kittel §5.2; Debye 1912]
        v_sound (m/s), n_density: atom density (m^-3). Returns K.
        """
        return (Constants.H_BAR / Constants.K_B) * v_sound * (6 * math.pi**2 * n_density)**(1/3)


# ─────────────────────────────────────────────────────────────────────────────
# UNIT HELPERS
# ─────────────────────────────────────────────────────────────────────────────

def km_to_m(km):          return km * 1e3
def m_to_km(m):           return m / 1e3
def au_to_m(au):          return au * Constants.AU
def m_to_au(m):           return m / Constants.AU
def ly_to_m(ly):          return ly * Constants.LY
def eV_to_J(eV):          return eV * Constants.E
def J_to_eV(J):           return J / Constants.E
def seconds_to_hours(s):  return s / 3600
def seconds_to_days(s):   return s / 86400
def degrees_to_rad(d):    return math.radians(d)
def rad_to_degrees(r):    return math.degrees(r)
def kg_to_solar(kg):      return kg / Constants.M_SUN
def solar_to_kg(ms):      return ms * Constants.M_SUN
def nm_to_m(nm):          return nm * 1e-9
def GHz_to_Hz(ghz):       return ghz * 1e9
def pc_to_m(pc):          return pc * Constants.PC
def m_to_pc(m):           return m / Constants.PC
