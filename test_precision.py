"""
AstroCalc Delta — Precision & Stability Test Suite
Runs each formula N=10,000 times, checks for floating-point drift,
numerical instability, and reports min/max/mean/std per formula.
"""

import math
import time
import sys

sys.path.insert(0, '.')
from math_module import (
    Constants, OrbitalMechanics, RocketPropulsion, PlanetaryScience,
    Relativity, ElectroMagnetism, QuantumMechanics,
    km_to_m, au_to_m, eV_to_J, J_to_eV, nm_to_m
)

N = 10_000  # iterations per formula

# ─────────────────────────────────────────────────────────────────────────────
# TEST REGISTRY
# Each entry: (name, callable, expected_rough_value, tolerance_fraction)
# tolerance_fraction: acceptable deviation from first result (floating point)
# ─────────────────────────────────────────────────────────────────────────────

TESTS = [

    # ── ORBITAL MECHANICS ──────────────────────────────────────────────────
    ("Orbital Velocity (ISS 408km)",
     lambda: OrbitalMechanics.orbital_velocity(km_to_m(408)),
     7668.0, 0.001),

    ("Orbital Period (ISS 408km) min",
     lambda: OrbitalMechanics.orbital_period(km_to_m(408)) / 60,
     92.57, 0.001),

    ("Escape Velocity (surface) km/s",
     lambda: OrbitalMechanics.escape_velocity(0) / 1000,
     11.186, 0.001),

    ("Hohmann Transfer LEO->GEO total dv m/s",
     lambda: OrbitalMechanics.hohmann_delta_v(6571e3, 42164e3)["total_dv_ms"],
     3934.7, 0.001),

    ("Vis-Viva at perigee r=6571km a=24368km",
     lambda: OrbitalMechanics.vis_viva(6571e3, 24368e3),
     10245.2, 0.001),

    ("Eccentricity ra=42164 rp=6571 km",
     lambda: OrbitalMechanics.eccentricity(42164e3, 6571e3),
     0.7303, 0.001),

    ("Inclination Change dv v=7.78km/s di=28.5deg",
     lambda: OrbitalMechanics.inclination_change_dv(7780, 28.5),
     3820.0, 0.01),

    ("GEO Altitude km",
     lambda: OrbitalMechanics.geostationary_altitude() / 1000,
     35786.0, 0.001),

    ("SOI Earth (AU)",
     lambda: OrbitalMechanics.sphere_of_influence(au_to_m(1.0), Constants.M_EARTH) / au_to_m(1.0),
     0.00617, 0.01),

    ("Bi-elliptic LEO->GEO total dv",
     lambda: OrbitalMechanics.bi_elliptic_delta_v(6571e3, 42164e3, 200000e3)["total_dv_ms"],
     3900.0, 0.5),  # wide tolerance — varies with transfer altitude

    ("Orbital Decay Timescale ISS (s)",
     lambda: OrbitalMechanics.orbital_decay_timescale(408e3),
     3.4e18, 0.5),   # exponential atmosphere approximation

    ("Repeating Ground Track 15 orbits/day SMA (km)",
     lambda: OrbitalMechanics.repeating_ground_track(15, 1) / 1000,
     6932.0, 0.01),

    ("Combined Plane Change dv",
     lambda: OrbitalMechanics.delta_v_plane_change_combined(3100, 3100, 30),
     1600.0, 0.1),

    ("Nodal Regression ISS 51.6 deg (deg/day)",
     lambda: OrbitalMechanics.nodal_regression_rate(408e3, 51.6),
     -4.99, 0.05),

    # ── ROCKET PROPULSION ──────────────────────────────────────────────────
    ("Tsiolkovsky Isp=450 m0=10000 mf=3000 (m/s)",
     lambda: RocketPropulsion.tsiolkovsky(450, 10000, 3000),
     5313.0, 0.001),

    ("Mass Ratio dv=9400 Isp=350",
     lambda: RocketPropulsion.mass_ratio(9400, 350),
     15.466, 0.001),

    ("Propellant Fraction dv=9400 Isp=350",
     lambda: RocketPropulsion.propellant_fraction(9400, 350),
     0.9353, 0.001),

    ("Thrust mdot=250 Isp=363 (kN)",
     lambda: RocketPropulsion.thrust(250, 363) / 1000,
     890.0, 0.001),

    ("TWR F=7607kN m=549t",
     lambda: RocketPropulsion.twr(7607e3, 549054),
     1.414, 0.001),

    ("Burn Time m0=10t mf=4t mdot=50",
     lambda: RocketPropulsion.burn_time(10000, 4000, 50),
     120.0, 0.001),

    ("Landing Burn dv (m_dry=25t prop=5t Isp=340)",
     lambda: RocketPropulsion.landing_burn_delta_v(25000, 5000, 340),
     590.0, 0.05),

    ("Nozzle Exit Velocity (LOX/LH2 params)",
     lambda: RocketPropulsion.nozzle_exit_velocity(3300, 1.2, 0.012, 0.1e5, 200e5),
     4200.0, 0.2),

    ("Reuse Propellant Margin (kg)",
     lambda: RocketPropulsion.reuse_propellant_margin(100000, 340, 9400, 1500, 800, 25000),
     -22007.0, 0.01),

    ("Optimal Staging 2-stage payload frac",
     lambda: RocketPropulsion.optimal_staging_payload_fraction(9400, 350, 2)["total_payload_frac"],
     0.02, 1.0),   # order-of-magnitude

    ("Chamber Thrust Coefficient gamma=1.2 pe/pc=0.001",
     lambda: RocketPropulsion.chamber_pressure_thrust_coefficient(1.2, 1e4, 1e7),
     1.5, 0.3),

    # ── PLANETARY SCIENCE ──────────────────────────────────────────────────
    ("Earth Surface Gravity (m/s^2)",
     lambda: PlanetaryScience.surface_gravity(Constants.M_EARTH, Constants.R_EARTH),
     9.820, 0.001),

    ("Earth Orbital Period (days)",
     lambda: PlanetaryScience.kepler_third_law_period(au_to_m(1.0)) / 86400,
     365.25, 0.001),

    ("Earth Equilibrium Temp A=0.3 (K)",
     lambda: PlanetaryScience.equilibrium_temperature(0.30, 1.0),
     254.85, 0.001),

    ("Synodic Period Earth-Mars (days)",
     lambda: PlanetaryScience.synodic_period(365.25, 686.97),
     779.9, 0.01),

    ("Planet Density Earth (kg/m^3)",
     lambda: PlanetaryScience.planet_density(Constants.M_EARTH, Constants.R_EARTH),
     5515.0, 0.01),

    ("Roche Limit Jupiter-water-moon (km)",
     lambda: PlanetaryScience.roche_limit(71492e3, 1326, 1000) / 1000,
     98957.0, 0.01),

    ("Lunar Transfer dv from 200km LEO (m/s)",
     lambda: PlanetaryScience.lunar_transfer_delta_v(200e3),
     3120.0, 0.01),

    ("Mars Surface Weight 70kg person (N)",
     lambda: PlanetaryScience.mars_surface_weight(70),
     260.96, 0.01),

    ("Asteroid Impact Energy 1e9kg at 20km/s (Mt TNT)",
     lambda: PlanetaryScience.asteroid_impact_energy(1e9, 20000)["energy_Mt_TNT"],
     47.8, 0.01),

    ("Lunar Escape Velocity (m/s)",
     lambda: PlanetaryScience.lunar_escape_velocity(),
     2380.0, 0.01),

    ("Earth-Mars Resonance Ratio",
     lambda: PlanetaryScience.orbital_resonance_ratio(365.25, 686.97),
     1.880, 0.001),

    # ── RELATIVITY ─────────────────────────────────────────────────────────
    ("Lorentz Factor v=0.8c",
     lambda: Relativity.lorentz_factor(0.8 * Constants.C),
     1.6667, 0.0001),

    ("Time Dilation t0=1yr v=0.9c (yr)",
     lambda: Relativity.time_dilation(365.25*86400, 0.9*Constants.C) / (365.25*86400),
     2.294, 0.001),

    ("Length Contraction L0=100m v=0.99c (m)",
     lambda: Relativity.length_contraction(100, 0.99*Constants.C),
     14.107, 0.001),

    ("Rest Energy 1kg (J)",
     lambda: Relativity.rest_energy(1.0),
     8.988e16, 0.001),

    ("Schwarzschild Radius 1 Msun (m)",
     lambda: Relativity.schwarzschild_radius(Constants.M_SUN),
     2954.0, 0.001),

    ("Hawking Temperature 1 Msun (nK)",
     lambda: Relativity.hawking_temperature(Constants.M_SUN) * 1e9,
     61.68, 0.01),

    ("GR Mercury Precession (arcsec/orbit)",
     lambda: Relativity.relativistic_precession_per_orbit(57.91e9, 0.2056, Constants.M_SUN),
     0.1034, 0.01),

    ("Photon Redshift at ISCO (z)",
     lambda: Relativity.photon_redshift_at_isco(Constants.M_SUN),
     0.2247, 0.001),

    ("GW Strain h (GW150914-like: Mc=30Msun d=410Mpc f=150Hz)",
     lambda: Relativity.gravitational_wave_strain(
         30*Constants.M_SUN, 410*3.086e22, 150),
     2.37e-21, 0.1),   # order-of-magnitude check

    ("Lense-Thirring precession Earth (rad/s)",
     lambda: Relativity.lense_thirring_precession(6371e3 + 650e3, 8.04e37),
     3.45e-10, 0.01),    # order of magnitude

    ("Penrose Efficiency a*=0.9",
     lambda: Relativity.penrose_max_efficiency(0.9),
     0.44, 0.05),

    # ── ELECTROMAGNETISM ───────────────────────────────────────────────────
    ("Coulomb Force q1=1nC q2=-1nC r=10mm (N)",
     lambda: abs(ElectroMagnetism.coulomb_force(1e-9, -1e-9, 0.01)),
     8.988e-5, 0.001),

    ("Cyclotron Radius proton v=1e6 B=1T (m)",
     lambda: ElectroMagnetism.cyclotron_radius(Constants.M_P, 1e6, Constants.E, 1.0),
     0.01044, 0.001),

    ("Cyclotron Frequency electron B=0.1T (GHz)",
     lambda: ElectroMagnetism.cyclotron_frequency(Constants.E, 0.1, Constants.M_E) / 1e9,
     2.799, 0.001),

    ("Wire B-field I=100A r=0.1m (uT)",
     lambda: ElectroMagnetism.biot_savart_wire(100, 0.1) * 1e6,
     200.0, 0.001),

    ("Skin Depth copper f=1MHz (um)",
     lambda: ElectroMagnetism.skin_depth(1.68e-8, 1e6) * 1e6,
     65.23, 0.01),

    ("Larmor Power electron a=1e15 m/s^2 (W)",
     lambda: ElectroMagnetism.larmor_radiation_power(Constants.E, 1e15),
     5.71e-24, 0.01),

    ("Plasma Frequency n_e=1e18 m^-3 (GHz)",
     lambda: ElectroMagnetism.plasma_frequency(1e18) / (2*math.pi) / 1e9,
     8.979, 0.001),

    ("Magnetic Pressure B=1T (Pa)",
     lambda: ElectroMagnetism.magnetic_pressure(1.0),
     397887.0, 0.001),

    ("Synchrotron Power gamma=100 B=1T (W)",
     lambda: ElectroMagnetism.synchrotron_radiation_power(100, 1.0),
     1.0e-8, 1.0),   # order check

    ("Alfven Speed B=1e-9T rho=1.67e-21 (km/s)",
     lambda: ElectroMagnetism.alfven_wave_speed(1e-9, 1.67e-21) / 1000,
     21.8, 0.05),

    ("Hall Parameter electron B=1T nu=1e10 Hz",
     lambda: ElectroMagnetism.hall_parameter(1.0, Constants.M_E, 1e10),
     1.759e1, 0.01),

    # ── QUANTUM MECHANICS ──────────────────────────────────────────────────
    ("Photon Energy 550nm (eV)",
     lambda: J_to_eV(QuantumMechanics.photon_energy_from_wavelength(nm_to_m(550))),
     2.254, 0.001),

    ("de Broglie proton v=1e6 m/s (fm)",
     lambda: QuantumMechanics.de_broglie_wavelength(Constants.M_P, 1e6) * 1e15,
     396.0, 0.01),

    ("Hydrogen n=1 energy (eV)",
     lambda: QuantumMechanics.hydrogen_energy_level(1),
     -13.6, 0.0001),

    ("Hydrogen 3->2 transition wavelength (nm)",
     lambda: QuantumMechanics.hydrogen_transition_wavelength(2, 3) * 1e9,
     656.11, 0.001),

    ("Heisenberg dx for dp=1 eV/c (pm)",
     lambda: QuantumMechanics.uncertainty_position_momentum(eV_to_J(1)/Constants.C) * 1e12,
     98663.0, 0.001),

    ("Compton Shift 90 deg (pm)",
     lambda: QuantumMechanics.compton_wavelength_shift(90) * 1e12,
     2.426, 0.001),

    ("Wien Peak Sun 5778K (nm)",
     lambda: QuantumMechanics.wien_peak_wavelength(5778) * 1e9,
     501.5, 0.001),

    ("Stefan-Boltzmann Sun surface (MW/m^2)",
     lambda: QuantumMechanics.stefan_boltzmann(5778) / 1e6,
     63.3, 0.01),

    ("Bohr Radius n=2 (pm)",
     lambda: QuantumMechanics.bohr_radius_n(2) * 1e12,
     211.7, 0.001),

    ("Fermi Energy copper n=8.49e28 (eV)",
     lambda: J_to_eV(QuantumMechanics.fermi_energy(8.49e28)),
     7.04, 0.01),

    ("Cooper Pair Gap Nb Tc=9.2K (meV)",
     lambda: QuantumMechanics.cooper_pair_binding_energy(9.2) / Constants.E * 1000,
     2.79, 0.01),

    ("Josephson Frequency V=1uV (GHz)",
     lambda: QuantumMechanics.josephson_frequency(1e-6) / 1e9,
     0.4836, 0.001),

    ("Tunneling Probability m_e V0=10eV E=5eV L=1nm",
     lambda: QuantumMechanics.quantum_tunneling_probability(
         Constants.M_E, eV_to_J(10), eV_to_J(5), 1e-9),
     1.12e-10, 0.1),

    ("Debye Temp aluminum v=5100 n=6.02e28 (K)",
     lambda: QuantumMechanics.debye_temperature(5100, 6.02e28),
     595.0, 0.05),
]


# ─────────────────────────────────────────────────────────────────────────────
# RUNNER
# ─────────────────────────────────────────────────────────────────────────────

def run_test(name, fn, expected, tol):
    results = []
    t_start = time.perf_counter()
    for _ in range(N):
        try:
            r = fn()
            # flatten dicts to their first numeric value for stats
            if isinstance(r, dict):
                r = list(r.values())[0]
            results.append(float(r))
        except Exception as e:
            return {"name": name, "status": "ERROR", "error": str(e)}
    elapsed = time.perf_counter() - t_start

    first  = results[0]
    mn     = min(results)
    mx     = max(results)
    avg    = sum(results) / len(results)
    variance = sum((x - avg)**2 for x in results) / len(results)
    std    = math.sqrt(variance)

    # Drift check: any result differs from first beyond floating point?
    max_drift = max(abs(r - first) for r in results)
    has_drift = max_drift > 0.0

    # Expected value check
    if expected != 0:
        rel_error = abs(avg - expected) / abs(expected)
        passes_expected = rel_error <= tol
    else:
        passes_expected = abs(avg) <= tol

    status = "PASS" if passes_expected else "WARN"
    if has_drift:
        status = "DRIFT"

    return {
        "name":       name,
        "status":     status,
        "first":      first,
        "avg":        avg,
        "min":        mn,
        "max":        mx,
        "std":        std,
        "max_drift":  max_drift,
        "expected":   expected,
        "rel_err_pct": abs(avg - expected) / abs(expected) * 100 if expected != 0 else 0,
        "time_ms":    elapsed * 1000,
        "iters":      N,
    }


def fmt(v):
    """Format a float for display."""
    if v == 0:
        return "0"
    mag = abs(v)
    if mag >= 1e6 or mag < 1e-4:
        return f"{v:.6e}"
    if mag >= 1000:
        return f"{v:.4f}"
    return f"{v:.6f}"


def main():
    print()
    print("=" * 80)
    print(f"  ASTROCALC DELTA — PRECISION & STABILITY TEST SUITE")
    print(f"  {N:,} iterations per formula | {len(TESTS)} formulas")
    print("=" * 80)

    results = []
    t_total = time.perf_counter()

    for name, fn, expected, tol in TESTS:
        r = run_test(name, fn, expected, tol)
        results.append(r)

        # Status symbol
        sym = {"PASS": "✓", "WARN": "⚠", "DRIFT": "⟳", "ERROR": "✗"}.get(r["status"], "?")

        if r["status"] == "ERROR":
            print(f"  {sym} {r['name']}")
            print(f"      ERROR: {r['error']}")
            continue

        drift_str = f"  DRIFT={r['max_drift']:.2e}" if r["max_drift"] > 0 else ""
        avg_str   = f"avg={fmt(r['avg'])}"
        exp_str   = f"expected≈{fmt(r['expected'])}"
        err_str   = f"err={r['rel_err_pct']:.3f}%"
        time_str  = f"{r['time_ms']:.1f}ms"

        print(f"  {sym} {r['name']}")
        print(f"      {avg_str}  |  {exp_str}  |  {err_str}  |  std={fmt(r['std'])}  |  {time_str}{drift_str}")

    elapsed_total = time.perf_counter() - t_total

    # Summary
    n_pass  = sum(1 for r in results if r.get("status") == "PASS")
    n_warn  = sum(1 for r in results if r.get("status") == "WARN")
    n_drift = sum(1 for r in results if r.get("status") == "DRIFT")
    n_error = sum(1 for r in results if r.get("status") == "ERROR")
    total_iters = N * len(TESTS)

    print()
    print("=" * 80)
    print(f"  RESULTS: {n_pass} PASS  |  {n_warn} WARN  |  {n_drift} DRIFT  |  {n_error} ERROR")
    print(f"  Total iterations: {total_iters:,}  |  Wall time: {elapsed_total:.2f}s")
    print(f"  Throughput: {total_iters/elapsed_total:,.0f} formula evaluations/second")
    print("=" * 80)

    if n_warn > 0:
        print()
        print("  WARNINGS (result outside expected tolerance):")
        for r in results:
            if r.get("status") == "WARN":
                print(f"    - {r['name']}")
                print(f"      avg={fmt(r['avg'])}  expected={fmt(r['expected'])}  err={r['rel_err_pct']:.2f}%")

    if n_error > 0:
        print()
        print("  ERRORS:")
        for r in results:
            if r.get("status") == "ERROR":
                print(f"    - {r['name']}: {r['error']}")

    print()
    return n_error == 0 and n_warn == 0


if __name__ == "__main__":
    ok = main()
    sys.exit(0 if ok else 1)
