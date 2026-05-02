"""
AstroCalc Delta — Flask API Server
Replaces the old Pyto/Django/HTTPServer setup.
Run locally:  python app.py
Deploy:       Railway / Render / Fly.io (Gunicorn via Procfile)
"""

import math
import os
import traceback
from flask import Flask, jsonify, request, send_from_directory
try:
    from flask_cors import CORS
    HAS_CORS = True
except ImportError:
    HAS_CORS = False

# Import all physics modules
from math_module import (
    Constants, OrbitalMechanics, RocketPropulsion,
    PlanetaryScience, Relativity, ElectroMagnetism, QuantumMechanics,
    km_to_m, au_to_m, eV_to_J, J_to_eV,
    seconds_to_hours, seconds_to_days, nm_to_m
)

app = Flask(__name__, static_folder="static", template_folder="templates")
if HAS_CORS:
    CORS(app)

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


# ─────────────────────────────────────────────────────────────────────────────
# SERVE FRONTEND
# ─────────────────────────────────────────────────────────────────────────────

@app.route("/")
def index():
    return send_from_directory(BASE_DIR, "index.html")

@app.route("/<path:filename>")
def static_files(filename):
    return send_from_directory(BASE_DIR, filename)


# ─────────────────────────────────────────────────────────────────────────────
# FORMULA REGISTRY
# ─────────────────────────────────────────────────────────────────────────────

FORMULA_REGISTRY = {

    # ── ORBITAL MECHANICS ────────────────────────────────────────────────────
    "orbital_velocity": {
        "label": "Orbital Velocity",
        "domain": "orbital",
        "equation": "v = √(μ / (R + h))",
        "params": [
            {"id": "altitude_km", "label": "Altitude", "unit": "km", "default": 400}
        ],
        "fn": lambda p: {
            "Orbital Velocity": f"{OrbitalMechanics.orbital_velocity(km_to_m(p['altitude_km'])):.2f} m/s",
            "Orbital Velocity (km/s)": f"{OrbitalMechanics.orbital_velocity(km_to_m(p['altitude_km']))/1000:.4f} km/s",
        }
    },
    "orbital_period": {
        "label": "Orbital Period",
        "domain": "orbital",
        "equation": "T = 2π√(r³/μ)",
        "params": [
            {"id": "altitude_km", "label": "Altitude", "unit": "km", "default": 400}
        ],
        "fn": lambda p: {
            "Period (s)": f"{OrbitalMechanics.orbital_period(km_to_m(p['altitude_km'])):.2f} s",
            "Period (min)": f"{OrbitalMechanics.orbital_period(km_to_m(p['altitude_km']))/60:.2f} min",
        }
    },
    "escape_velocity": {
        "label": "Escape Velocity",
        "domain": "orbital",
        "equation": "v_esc = √(2μ / r)",
        "params": [
            {"id": "altitude_km", "label": "Altitude", "unit": "km", "default": 0}
        ],
        "fn": lambda p: {
            "Escape Velocity": f"{OrbitalMechanics.escape_velocity(km_to_m(p['altitude_km'])):.2f} m/s",
            "Escape Velocity (km/s)": f"{OrbitalMechanics.escape_velocity(km_to_m(p['altitude_km']))/1000:.3f} km/s",
        }
    },
    "hohmann_transfer": {
        "label": "Hohmann Transfer",
        "domain": "orbital",
        "equation": "Δv_total = Δv₁ + Δv₂",
        "params": [
            {"id": "r1_km", "label": "Orbit 1 Radius", "unit": "km", "default": 6771},
            {"id": "r2_km", "label": "Orbit 2 Radius", "unit": "km", "default": 42164}
        ],
        "fn": lambda p: {k: v for k, v in {
            **{k2: f"{v2} m/s" for k2, v2 in OrbitalMechanics.hohmann_delta_v(
                km_to_m(p['r1_km']), km_to_m(p['r2_km'])).items()}
        }.items()}
    },
    "vis_viva": {
        "label": "Vis-Viva Velocity",
        "domain": "orbital",
        "equation": "v = √(μ(2/r − 1/a))",
        "params": [
            {"id": "r_km", "label": "Current Radius", "unit": "km", "default": 6771},
            {"id": "a_km", "label": "Semi-major Axis", "unit": "km", "default": 24468}
        ],
        "fn": lambda p: {
            "Velocity": f"{OrbitalMechanics.vis_viva(km_to_m(p['r_km']), km_to_m(p['a_km'])):.2f} m/s"
        }
    },
    "eccentricity": {
        "label": "Orbital Eccentricity",
        "domain": "orbital",
        "equation": "e = (r_a − r_p) / (r_a + r_p)",
        "params": [
            {"id": "r_apo_km", "label": "Apoapsis Radius", "unit": "km", "default": 42164},
            {"id": "r_per_km", "label": "Periapsis Radius", "unit": "km", "default": 6571}
        ],
        "fn": lambda p: {
            "Eccentricity": f"{OrbitalMechanics.eccentricity(km_to_m(p['r_apo_km']), km_to_m(p['r_per_km'])):.6f}"
        }
    },
    "inclination_change": {
        "label": "Inclination Change Δv",
        "domain": "orbital",
        "equation": "Δv = 2v·sin(Δi/2)",
        "params": [
            {"id": "v_kms", "label": "Orbital Velocity", "unit": "km/s", "default": 7.78},
            {"id": "delta_i_deg", "label": "Inclination Change", "unit": "°", "default": 28.5}
        ],
        "fn": lambda p: {
            "Δv Required": f"{OrbitalMechanics.inclination_change_dv(p['v_kms']*1000, p['delta_i_deg']):.2f} m/s"
        }
    },
    "geostationary_altitude": {
        "label": "GEO Altitude",
        "domain": "orbital",
        "equation": "r_GEO = (μ/ω²)^(1/3) − R",
        "params": [],
        "fn": lambda p: {
            "GEO Altitude": f"{OrbitalMechanics.geostationary_altitude()/1000:.2f} km",
            "GEO Radius": f"{(OrbitalMechanics.geostationary_altitude()+Constants.R_EARTH)/1000:.2f} km",
        }
    },

    # ── ROCKET PROPULSION ────────────────────────────────────────────────────
    "tsiolkovsky": {
        "label": "Tsiolkovsky Rocket Equation",
        "domain": "propulsion",
        "equation": "Δv = Isp·g₀·ln(m₀/mf)",
        "params": [
            {"id": "isp", "label": "Specific Impulse (Isp)", "unit": "s", "default": 450},
            {"id": "m0_kg", "label": "Initial Mass (m₀)", "unit": "kg", "default": 10000},
            {"id": "mf_kg", "label": "Final Mass (mf)", "unit": "kg", "default": 3000}
        ],
        "fn": lambda p: {
            "Delta-v": f"{RocketPropulsion.tsiolkovsky(p['isp'], p['m0_kg'], p['mf_kg']):.2f} m/s",
            "Delta-v (km/s)": f"{RocketPropulsion.tsiolkovsky(p['isp'], p['m0_kg'], p['mf_kg'])/1000:.4f} km/s",
        }
    },
    "mass_ratio": {
        "label": "Required Mass Ratio",
        "domain": "propulsion",
        "equation": "m₀/mf = e^(Δv/ve)",
        "params": [
            {"id": "delta_v_kms", "label": "Delta-v", "unit": "km/s", "default": 9.4},
            {"id": "isp", "label": "Isp", "unit": "s", "default": 350}
        ],
        "fn": lambda p: {
            "Mass Ratio (m₀/mf)": f"{RocketPropulsion.mass_ratio(p['delta_v_kms']*1000, p['isp']):.4f}",
            "Propellant Fraction": f"{RocketPropulsion.propellant_fraction(p['delta_v_kms']*1000, p['isp'])*100:.2f} %"
        }
    },
    "thrust": {
        "label": "Engine Thrust",
        "domain": "propulsion",
        "equation": "F = ṁ·Isp·g₀",
        "params": [
            {"id": "mdot_kgs", "label": "Mass Flow Rate", "unit": "kg/s", "default": 250},
            {"id": "isp", "label": "Isp", "unit": "s", "default": 363}
        ],
        "fn": lambda p: {
            "Thrust": f"{RocketPropulsion.thrust(p['mdot_kgs'], p['isp'])/1000:.2f} kN",
            "Exhaust Velocity": f"{RocketPropulsion.exhaust_velocity(p['isp']):.2f} m/s"
        }
    },
    "twr": {
        "label": "Thrust-to-Weight Ratio",
        "domain": "propulsion",
        "equation": "TWR = F / (m·g₀)",
        "params": [
            {"id": "thrust_kn", "label": "Thrust", "unit": "kN", "default": 7607},
            {"id": "mass_kg", "label": "Total Mass", "unit": "kg", "default": 549054}
        ],
        "fn": lambda p: {
            "TWR": f"{RocketPropulsion.twr(p['thrust_kn']*1000, p['mass_kg']):.4f}"
        }
    },
    "burn_time": {
        "label": "Engine Burn Time",
        "domain": "propulsion",
        "equation": "t = (m₀ − mf) / ṁ",
        "params": [
            {"id": "m0_kg", "label": "Initial Mass", "unit": "kg", "default": 10000},
            {"id": "mf_kg", "label": "Final Mass", "unit": "kg", "default": 4000},
            {"id": "mdot_kgs", "label": "Mass Flow Rate", "unit": "kg/s", "default": 50}
        ],
        "fn": lambda p: {
            "Burn Time": f"{RocketPropulsion.burn_time(p['m0_kg'], p['mf_kg'], p['mdot_kgs']):.2f} s",
            "Burn Time (min)": f"{RocketPropulsion.burn_time(p['m0_kg'], p['mf_kg'], p['mdot_kgs'])/60:.2f} min"
        }
    },

    "landing_burn_dv": {
        "label": "Landing Burn Δv",
        "domain": "propulsion",
        "equation": "Δv = Isp·g₀·ln((m_dry+m_prop)/m_dry)",
        "params": [
            {"id": "mass_dry_kg", "label": "Dry Mass", "unit": "kg", "default": 25000},
            {"id": "mass_prop_kg", "label": "Landing Propellant", "unit": "kg", "default": 5000},
            {"id": "isp", "label": "Isp", "unit": "s", "default": 340}
        ],
        "fn": lambda p: {
            "Landing Burn Δv": f"{RocketPropulsion.landing_burn_delta_v(p['mass_dry_kg'], p['mass_prop_kg'], p['isp']):.2f} m/s",
            "Landing Burn Δv (km/s)": f"{RocketPropulsion.landing_burn_delta_v(p['mass_dry_kg'], p['mass_prop_kg'], p['isp'])/1000:.4f} km/s"
        }
    },
    "nozzle_exit_velocity": {
        "label": "Nozzle Exit Velocity",
        "domain": "propulsion",
        "equation": "ve = √(2γ/(γ−1)·RT_c/M·[1−(pe/pc)^((γ−1)/γ)])",
        "params": [
            {"id": "T_chamber_K", "label": "Chamber Temperature", "unit": "K", "default": 3300},
            {"id": "gamma", "label": "Ratio of Specific Heats γ", "unit": "", "default": 1.2},
            {"id": "M_mol", "label": "Molar Mass", "unit": "kg/mol", "default": 0.012},
            {"id": "p_e_bar", "label": "Exit Pressure", "unit": "bar", "default": 0.1},
            {"id": "p_c_bar", "label": "Chamber Pressure", "unit": "bar", "default": 200}
        ],
        "fn": lambda p: {
            "Exit Velocity": f"{RocketPropulsion.nozzle_exit_velocity(p['T_chamber_K'], p['gamma'], p['M_mol'], p['p_e_bar']*1e5, p['p_c_bar']*1e5):.2f} m/s",
            "Effective Isp": f"{RocketPropulsion.nozzle_exit_velocity(p['T_chamber_K'], p['gamma'], p['M_mol'], p['p_e_bar']*1e5, p['p_c_bar']*1e5)/Constants.G_0:.1f} s"
        }
    },
    "reuse_margin": {
        "label": "Reuse Propellant Margin",
        "domain": "propulsion",
        "equation": "m_margin = m₀·e^(−Δv_total/ve) − m_dry",
        "params": [
            {"id": "m0_kg", "label": "Liftoff Mass", "unit": "kg", "default": 549054},
            {"id": "isp", "label": "Isp", "unit": "s", "default": 340},
            {"id": "dv_ascent_kms", "label": "Ascent Δv", "unit": "km/s", "default": 2.0},
            {"id": "dv_reentry_kms", "label": "Reentry Burn Δv", "unit": "km/s", "default": 0.4},
            {"id": "dv_landing_kms", "label": "Landing Burn Δv", "unit": "km/s", "default": 0.25},
            {"id": "mf_dry_kg", "label": "Dry Mass", "unit": "kg", "default": 25000}
        ],
        "fn": lambda p: {
            "Propellant Margin": f"{RocketPropulsion.reuse_propellant_margin(p['m0_kg'], p['isp'], p['dv_ascent_kms']*1000, p['dv_reentry_kms']*1000, p['dv_landing_kms']*1000, p['mf_dry_kg']):.2f} kg",
            "Reusable?": "Yes ✓" if RocketPropulsion.reuse_propellant_margin(p['m0_kg'], p['isp'], p['dv_ascent_kms']*1000, p['dv_reentry_kms']*1000, p['dv_landing_kms']*1000, p['mf_dry_kg']) > 0 else "No — insufficient prop"
        }
    },
    "optimal_staging": {
        "label": "Optimal Staging",
        "domain": "propulsion",
        "equation": "PF = ((1−ε)/MR − ε)ⁿ",
        "params": [
            {"id": "dv_total_kms", "label": "Total Δv", "unit": "km/s", "default": 9.4},
            {"id": "isp", "label": "Isp", "unit": "s", "default": 350},
            {"id": "n_stages", "label": "Number of Stages", "unit": "", "default": 2},
            {"id": "eps", "label": "Structural Fraction ε", "unit": "", "default": 0.1}
        ],
        "fn": lambda p: (lambda r: {
            "Mass Ratio per Stage": f"{r['mass_ratio_per_stage']:.4f}",
            "Payload Frac per Stage": f"{r['payload_frac_per_stage']*100:.4f} %",
            "Total Payload Fraction": f"{r['total_payload_frac']*100:.6f} %"
        })(RocketPropulsion.optimal_staging_payload_fraction(p['dv_total_kms']*1000, p['isp'], int(p['n_stages']), p['eps']))
    },
    "gravity_loss": {
        "label": "Gravity Loss Estimate",
        "domain": "propulsion",
        "equation": "Δv_loss ≈ g₀·t·(1 − 1/TWR)",
        "params": [
            {"id": "twr", "label": "Thrust-to-Weight Ratio", "unit": "", "default": 1.5},
            {"id": "burn_time_s", "label": "Burn Duration", "unit": "s", "default": 150}
        ],
        "fn": lambda p: {
            "Gravity Loss": f"{RocketPropulsion.gravity_loss_estimate(p['twr'], p['burn_time_s']):.2f} m/s",
            "Gravity Loss (km/s)": f"{RocketPropulsion.gravity_loss_estimate(p['twr'], p['burn_time_s'])/1000:.4f} km/s"
        }
    },

    # ── PLANETARY SCIENCE ────────────────────────────────────────────────────
    "surface_gravity": {
        "label": "Surface Gravity",
        "domain": "planetary",
        "equation": "g = GM / R²",
        "params": [
            {"id": "mass_kg_exp", "label": "Body Mass (×10²⁴ kg)", "unit": "×10²⁴ kg", "default": 5.972},
            {"id": "radius_km", "label": "Body Radius", "unit": "km", "default": 6371}
        ],
        "fn": lambda p: {
            "Surface Gravity": f"{PlanetaryScience.surface_gravity(p['mass_kg_exp']*1e24, km_to_m(p['radius_km'])):.4f} m/s²",
            "vs Earth g": f"{PlanetaryScience.surface_gravity(p['mass_kg_exp']*1e24, km_to_m(p['radius_km']))/9.80665:.3f} g"
        }
    },
    "kepler_period": {
        "label": "Orbital Period (Kepler III)",
        "domain": "planetary",
        "equation": "T = 2π√(a³/μ_sun)",
        "params": [
            {"id": "a_au", "label": "Semi-major Axis", "unit": "AU", "default": 1.0}
        ],
        "fn": lambda p: {
            "Period (s)": f"{PlanetaryScience.kepler_third_law_period(au_to_m(p['a_au'])):.4e} s",
            "Period (days)": f"{seconds_to_days(PlanetaryScience.kepler_third_law_period(au_to_m(p['a_au']))):.4f} days",
            "Period (years)": f"{seconds_to_days(PlanetaryScience.kepler_third_law_period(au_to_m(p['a_au'])))/365.25:.4f} yr"
        }
    },
    "roche_limit": {
        "label": "Roche Limit",
        "domain": "planetary",
        "equation": "d = R_p·(2ρ_p/ρ_s)^(1/3)",
        "params": [
            {"id": "R_primary_km", "label": "Primary Radius", "unit": "km", "default": 71492},
            {"id": "rho_primary", "label": "Primary Density", "unit": "kg/m³", "default": 1326},
            {"id": "rho_satellite", "label": "Satellite Density", "unit": "kg/m³", "default": 1000}
        ],
        "fn": lambda p: {
            "Roche Limit": f"{PlanetaryScience.roche_limit(km_to_m(p['R_primary_km']), p['rho_primary'], p['rho_satellite'])/1000:.2f} km"
        }
    },
    "equilibrium_temp": {
        "label": "Planetary Equilibrium Temperature",
        "domain": "planetary",
        "equation": "T = T_sun·√(R_sun/2d)·(1−A)^¼",
        "params": [
            {"id": "albedo", "label": "Bond Albedo", "unit": "", "default": 0.30},
            {"id": "distance_au", "label": "Distance from Star", "unit": "AU", "default": 1.0}
        ],
        "fn": lambda p: {
            "Equilibrium Temp": f"{PlanetaryScience.equilibrium_temperature(p['albedo'], p['distance_au']):.2f} K",
            "Equilibrium Temp (°C)": f"{PlanetaryScience.equilibrium_temperature(p['albedo'], p['distance_au'])-273.15:.2f} °C"
        }
    },
    "synodic_period": {
        "label": "Synodic Period",
        "domain": "planetary",
        "equation": "1/T_syn = |1/T₁ − 1/T₂|",
        "params": [
            {"id": "T1_days", "label": "Inner Period", "unit": "days", "default": 365.25},
            {"id": "T2_days", "label": "Outer Period", "unit": "days", "default": 686.97}
        ],
        "fn": lambda p: {
            "Synodic Period": f"{PlanetaryScience.synodic_period(p['T1_days'], p['T2_days']):.4f} days"
        }
    },
    "planet_density": {
        "label": "Mean Planetary Density",
        "domain": "planetary",
        "equation": "ρ = M / (4/3 π R³)",
        "params": [
            {"id": "mass_kg_exp", "label": "Mass (×10²⁴ kg)", "unit": "×10²⁴ kg", "default": 5.972},
            {"id": "radius_km", "label": "Radius", "unit": "km", "default": 6371}
        ],
        "fn": lambda p: {
            "Mean Density": f"{PlanetaryScience.planet_density(p['mass_kg_exp']*1e24, km_to_m(p['radius_km'])):.2f} kg/m³"
        }
    },

    # ── RELATIVITY & BLACK HOLES ─────────────────────────────────────────────
    "lorentz_factor": {
        "label": "Lorentz Factor",
        "domain": "relativity",
        "equation": "γ = 1 / √(1 − v²/c²)",
        "params": [
            {"id": "v_frac_c", "label": "Velocity (fraction of c)", "unit": "c", "default": 0.5}
        ],
        "fn": lambda p: {
            "γ (Lorentz Factor)": f"{Relativity.lorentz_factor(p['v_frac_c']*Constants.C):.6f}"
        }
    },
    "time_dilation": {
        "label": "Time Dilation",
        "domain": "relativity",
        "equation": "t = γ·t₀",
        "params": [
            {"id": "t0_s", "label": "Proper Time (t₀)", "unit": "s", "default": 1.0},
            {"id": "v_frac_c", "label": "Velocity (fraction of c)", "unit": "c", "default": 0.8}
        ],
        "fn": lambda p: {
            "Dilated Time": f"{Relativity.time_dilation(p['t0_s'], p['v_frac_c']*Constants.C):.6f} s",
            "Time Stretched By": f"{Relativity.lorentz_factor(p['v_frac_c']*Constants.C):.4f}×"
        }
    },
    "length_contraction": {
        "label": "Length Contraction",
        "domain": "relativity",
        "equation": "L = L₀ / γ",
        "params": [
            {"id": "L0_m", "label": "Proper Length (L₀)", "unit": "m", "default": 100.0},
            {"id": "v_frac_c", "label": "Velocity (fraction of c)", "unit": "c", "default": 0.9}
        ],
        "fn": lambda p: {
            "Contracted Length": f"{Relativity.length_contraction(p['L0_m'], p['v_frac_c']*Constants.C):.6f} m"
        }
    },
    "rest_energy": {
        "label": "Rest Energy (E=mc²)",
        "domain": "relativity",
        "equation": "E₀ = mc²",
        "params": [
            {"id": "mass_kg", "label": "Mass", "unit": "kg", "default": 1.0}
        ],
        "fn": lambda p: {
            "Rest Energy": f"{Relativity.rest_energy(p['mass_kg']):.4e} J",
            "Rest Energy (MeV)": f"{J_to_eV(Relativity.rest_energy(p['mass_kg']))/1e6:.4e} MeV"
        }
    },
    "schwarzschild_radius": {
        "label": "Schwarzschild Radius",
        "domain": "relativity",
        "equation": "rs = 2GM / c²",
        "params": [
            {"id": "mass_solar", "label": "Mass", "unit": "M_sun", "default": 1.0}
        ],
        "fn": lambda p: {
            "Schwarzschild Radius": f"{Relativity.schwarzschild_radius(p['mass_solar']*Constants.M_SUN):.4f} m",
            "Schwarzschild Radius (km)": f"{Relativity.schwarzschild_radius(p['mass_solar']*Constants.M_SUN)/1000:.4f} km"
        }
    },
    "hawking_temperature": {
        "label": "Hawking Temperature",
        "domain": "relativity",
        "equation": "T_H = ℏc³ / (8πGMk_B)",
        "params": [
            {"id": "mass_solar", "label": "Black Hole Mass", "unit": "M_sun", "default": 1.0}
        ],
        "fn": lambda p: {
            "Hawking Temperature": f"{Relativity.hawking_temperature(p['mass_solar']*Constants.M_SUN):.4e} K",
            "Evaporation Time": f"{Relativity.bh_evaporation_time(p['mass_solar']*Constants.M_SUN):.4e} s"
        }
    },
    "kerr_horizon": {
        "label": "Kerr Event Horizon",
        "domain": "relativity",
        "equation": "r₊ = GM/c² + √((GM/c²)² − a²)",
        "params": [
            {"id": "mass_solar", "label": "Black Hole Mass", "unit": "M_sun", "default": 10.0},
            {"id": "spin_param_frac", "label": "Spin Parameter a/a_max", "unit": "", "default": 0.5}
        ],
        "fn": lambda p: (lambda M, a_max: {
            "Outer Horizon r₊": f"{Relativity.kerr_event_horizon(M, p['spin_param_frac']*a_max):.4f} m" if Relativity.kerr_event_horizon(M, p['spin_param_frac']*a_max) else "Naked singularity",
            "Schwarzschild r_s": f"{Relativity.schwarzschild_radius(M):.4f} m",
            "Photon Sphere": f"{Relativity.photon_sphere_radius(M):.4f} m",
            "ISCO": f"{Relativity.innermost_stable_orbit(M):.4f} m"
        })(p['mass_solar']*Constants.M_SUN, Constants.G*p['mass_solar']*Constants.M_SUN/Constants.C**2)
    },
    "gravitational_redshift": {
        "label": "Gravitational Redshift",
        "domain": "relativity",
        "equation": "z+1 = 1/√(1 − rs/r)",
        "params": [
            {"id": "mass_solar", "label": "Body Mass", "unit": "M_sun", "default": 1.0},
            {"id": "r_rs_frac", "label": "Distance (multiples of rs)", "unit": "× rs", "default": 10.0}
        ],
        "fn": lambda p: (lambda rs: {
            "Redshift z": f"{Relativity.gravitational_redshift(p['r_rs_frac']*rs, p['mass_solar']*Constants.M_SUN)-1:.6f}",
            "Factor (z+1)": f"{Relativity.gravitational_redshift(p['r_rs_frac']*rs, p['mass_solar']*Constants.M_SUN):.6f}"
        })(Relativity.schwarzschild_radius(p['mass_solar']*Constants.M_SUN))
    },

    # ── ELECTROMAGNETISM ─────────────────────────────────────────────────────
    "coulomb_force": {
        "label": "Coulomb's Law",
        "domain": "em",
        "equation": "F = kq₁q₂ / r²",
        "params": [
            {"id": "q1_nC", "label": "Charge q₁", "unit": "nC", "default": 1.0},
            {"id": "q2_nC", "label": "Charge q₂", "unit": "nC", "default": -1.0},
            {"id": "r_mm", "label": "Separation r", "unit": "mm", "default": 10.0}
        ],
        "fn": lambda p: {
            "Force": f"{ElectroMagnetism.coulomb_force(p['q1_nC']*1e-9, p['q2_nC']*1e-9, p['r_mm']*1e-3):.6e} N"
        }
    },
    "cyclotron_radius": {
        "label": "Cyclotron Radius",
        "domain": "em",
        "equation": "r = mv / (qB)",
        "params": [
            {"id": "v_kms", "label": "Particle Velocity", "unit": "km/s", "default": 1000.0},
            {"id": "B_T", "label": "Magnetic Field B", "unit": "T", "default": 1e-5},
            {"id": "is_proton", "label": "Particle (1=proton, 0=electron)", "unit": "", "default": 1}
        ],
        "fn": lambda p: (lambda m: {
            "Cyclotron Radius": f"{ElectroMagnetism.cyclotron_radius(m, p['v_kms']*1000, Constants.E, p['B_T']):.4f} m",
            "Cyclotron Freq": f"{ElectroMagnetism.cyclotron_frequency(Constants.E, p['B_T'], m):.4f} Hz"
        })(Constants.M_P if p['is_proton'] >= 0.5 else Constants.M_E)
    },
    "larmor_power": {
        "label": "Larmor Radiation Power",
        "domain": "em",
        "equation": "P = q²a² / (6πε₀c³)",
        "params": [
            {"id": "acceleration_ms2", "label": "Acceleration", "unit": "m/s²", "default": 1e15}
        ],
        "fn": lambda p: {
            "Radiated Power (electron)": f"{ElectroMagnetism.larmor_radiation_power(Constants.E, p['acceleration_ms2']):.4e} W"
        }
    },

    # ── QUANTUM MECHANICS ────────────────────────────────────────────────────
    "photon_energy": {
        "label": "Photon Energy",
        "domain": "quantum",
        "equation": "E = hc / λ",
        "params": [
            {"id": "wavelength_nm", "label": "Wavelength", "unit": "nm", "default": 550.0}
        ],
        "fn": lambda p: {
            "Photon Energy (J)": f"{QuantumMechanics.photon_energy_from_wavelength(nm_to_m(p['wavelength_nm'])):.4e} J",
            "Photon Energy (eV)": f"{J_to_eV(QuantumMechanics.photon_energy_from_wavelength(nm_to_m(p['wavelength_nm']))):.4f} eV",
            "Frequency": f"{QuantumMechanics.photon_frequency(nm_to_m(p['wavelength_nm'])):.4e} Hz"
        }
    },
    "de_broglie": {
        "label": "de Broglie Wavelength",
        "domain": "quantum",
        "equation": "λ = h / mv",
        "params": [
            {"id": "mass_amu", "label": "Particle Mass", "unit": "amu", "default": 1.0},
            {"id": "v_kms", "label": "Velocity", "unit": "km/s", "default": 1.0}
        ],
        "fn": lambda p: {
            "de Broglie λ": f"{QuantumMechanics.de_broglie_wavelength(p['mass_amu']*1.66053906660e-27, p['v_kms']*1000):.4e} m",
            "de Broglie λ (pm)": f"{QuantumMechanics.de_broglie_wavelength(p['mass_amu']*1.66053906660e-27, p['v_kms']*1000)*1e12:.4f} pm"
        }
    },
    "hydrogen_levels": {
        "label": "Hydrogen Energy Levels",
        "domain": "quantum",
        "equation": "E_n = −13.6 eV / n²",
        "params": [
            {"id": "n_upper", "label": "Upper Level n", "unit": "", "default": 3},
            {"id": "n_lower", "label": "Lower Level m", "unit": "", "default": 2}
        ],
        "fn": lambda p: (lambda n1, n2: {
            "Upper Level Energy": f"{QuantumMechanics.hydrogen_energy_level(int(n1)):.4f} eV",
            "Lower Level Energy": f"{QuantumMechanics.hydrogen_energy_level(int(n2)):.4f} eV",
            "Photon Wavelength": f"{(QuantumMechanics.hydrogen_transition_wavelength(int(n2), int(n1)) or 0)*1e9:.4f} nm" if n1 > n2 else "N/A"
        })(p['n_upper'], p['n_lower'])
    },
    "heisenberg": {
        "label": "Heisenberg Uncertainty",
        "domain": "quantum",
        "equation": "Δx·Δp ≥ ℏ/2",
        "params": [
            {"id": "delta_p_eVc", "label": "Momentum Uncertainty (eV/c)", "unit": "eV/c", "default": 1.0}
        ],
        "fn": lambda p: {
            "Min Δx": f"{QuantumMechanics.uncertainty_position_momentum(eV_to_J(p['delta_p_eVc'])/Constants.C):.4e} m",
            "Min Δx (pm)": f"{QuantumMechanics.uncertainty_position_momentum(eV_to_J(p['delta_p_eVc'])/Constants.C)*1e12:.4f} pm"
        }
    },
    "wien_law": {
        "label": "Wien's Displacement Law",
        "domain": "quantum",
        "equation": "λ_max = b / T",
        "params": [
            {"id": "temp_K", "label": "Temperature", "unit": "K", "default": 5778}
        ],
        "fn": lambda p: {
            "Peak Wavelength": f"{QuantumMechanics.wien_peak_wavelength(p['temp_K'])*1e9:.2f} nm",
            "Radiated Power": f"{QuantumMechanics.stefan_boltzmann(p['temp_K']):.4e} W/m²"
        }
    },
    "compton_shift": {
        "label": "Compton Scattering",
        "domain": "quantum",
        "equation": "Δλ = (h/m_e·c)·(1−cosθ)",
        "params": [
            {"id": "theta_deg", "label": "Scattering Angle θ", "unit": "°", "default": 90.0}
        ],
        "fn": lambda p: {
            "Wavelength Shift Δλ": f"{QuantumMechanics.compton_wavelength_shift(p['theta_deg'])*1e12:.6f} pm",
            "Compton Wavelength λ_c": f"{Constants.H/(Constants.M_E*Constants.C)*1e12:.6f} pm"
        }
    },
}


# ─────────────────────────────────────────────────────────────────────────────
# API ROUTES
# ─────────────────────────────────────────────────────────────────────────────

@app.route("/api/formulas")
def get_formulas():
    """Return the full formula catalogue (no functions)."""
    catalogue = {}
    for key, f in FORMULA_REGISTRY.items():
        catalogue[key] = {
            "label": f["label"],
            "domain": f["domain"],
            "equation": f["equation"],
            "params": f["params"]
        }
    return jsonify(catalogue)


PRECISION_ITERS = 10_000

@app.route("/api/calculate", methods=["POST"])
def calculate():
    """Execute a formula N times, return results + precision stats."""
    import time
    data = request.get_json(force=True)
    formula_key = data.get("formula")
    raw_params  = data.get("params", {})

    if formula_key not in FORMULA_REGISTRY:
        return jsonify({"error": f"Unknown formula: {formula_key}"}), 404

    entry = FORMULA_REGISTRY[formula_key]

    try:
        params = {k: float(v) for k, v in raw_params.items()}

        # Run once to get formatted results
        result = entry["fn"](params)

        # Precision loop — extract first numeric value from result dict
        # and run N times to check for drift
        first_vals = {}
        all_nums   = {}
        for k, v in result.items():
            # Strip units — grab leading number
            try:
                num = float(str(v).split()[0].replace(',',''))
                first_vals[k] = num
                all_nums[k]   = [num]
            except Exception:
                pass

        t0 = time.perf_counter()
        for _ in range(PRECISION_ITERS - 1):
            r2 = entry["fn"](params)
            for k in all_nums:
                try:
                    num = float(str(r2[k]).split()[0].replace(',',''))
                    all_nums[k].append(num)
                except Exception:
                    pass
        elapsed_ms = (time.perf_counter() - t0) * 1000

        # Compute stats per output field
        precision = {}
        any_drift = False
        for k, vals in all_nums.items():
            if not vals:
                continue
            mn  = min(vals)
            mx  = max(vals)
            avg = sum(vals) / len(vals)
            drift = mx - mn
            if drift > 0:
                any_drift = True
            precision[k] = {
                "min": mn, "max": mx, "avg": avg,
                "drift": drift, "stable": drift == 0.0
            }

        return jsonify({
            "ok":            True,
            "results":       result,
            "iterations":    PRECISION_ITERS,
            "elapsed_ms":    round(elapsed_ms, 2),
            "any_drift":     any_drift,
            "precision":     precision,
        })

    except Exception as exc:
        traceback.print_exc()
        return jsonify({"ok": False, "error": str(exc)}), 400


@app.route("/api/constants")
def get_constants():
    """Return key physical constants."""
    return jsonify({
        "Speed of Light":       f"{Constants.C:.6e} m/s",
        "Gravitational Const":  f"{Constants.G:.6e} m³ kg⁻¹ s⁻²",
        "Planck Constant":      f"{Constants.H:.6e} J·s",
        "Boltzmann Constant":   f"{Constants.K_B:.6e} J/K",
        "Elementary Charge":    f"{Constants.E:.6e} C",
        "Electron Mass":        f"{Constants.M_E:.6e} kg",
        "Proton Mass":          f"{Constants.M_P:.6e} kg",
        "Earth Mass":           f"{Constants.M_EARTH:.6e} kg",
        "Earth Radius":         f"{Constants.R_EARTH:.6e} m",
        "Solar Mass":           f"{Constants.M_SUN:.6e} kg",
        "Solar Luminosity":     f"{Constants.L_SUN:.6e} W",
        "Astronomical Unit":    f"{Constants.AU:.6e} m",
        "Stefan-Boltzmann σ":   f"{Constants.SIGMA:.6e} W m⁻² K⁻⁴",
    })


@app.route("/api/health")
def health():
    return jsonify({"status": "AstroCalc Delta online"})


# ─────────────────────────────────────────────────────────────────────────────
# PLOT REGISTRY
# Each entry defines how to sweep a parameter and what curves to return.
# ─────────────────────────────────────────────────────────────────────────────

def linspace(start, stop, n=200):
    step = (stop - start) / (n - 1)
    return [start + i * step for i in range(n)]

PLOT_REGISTRY = {

    # ── ORBITAL ──────────────────────────────────────────────────────────────
    "orbital_velocity": {
        "title": "Orbital Velocity vs Altitude",
        "x_label": "Altitude (km)",
        "x_unit": "km",
        "curves": [
            {
                "label": "Orbital Velocity (km/s)",
                "color": "#3a7fff",
                "fn": lambda alt_km: OrbitalMechanics.orbital_velocity(km_to_m(alt_km)) / 1000
            },
            {
                "label": "Escape Velocity (km/s)",
                "color": "#ff4a6e",
                "fn": lambda alt_km: OrbitalMechanics.escape_velocity(km_to_m(alt_km)) / 1000
            }
        ],
        "x_range": [200, 36000],
        "annotations": [
            {"x": 408,   "label": "ISS (408 km)"},
            {"x": 35786, "label": "GEO (35,786 km)"},
        ]
    },

    "orbital_period": {
        "title": "Orbital Period vs Altitude",
        "x_label": "Altitude (km)",
        "x_unit": "km",
        "curves": [
            {
                "label": "Period (minutes)",
                "color": "#2affa0",
                "fn": lambda alt_km: OrbitalMechanics.orbital_period(km_to_m(alt_km)) / 60
            }
        ],
        "x_range": [200, 36000],
        "annotations": [
            {"x": 408,   "label": "ISS"},
            {"x": 35786, "label": "GEO"},
        ]
    },

    "hohmann_transfer": {
        "title": "Hohmann Transfer Δv vs Target Orbit",
        "x_label": "Target Orbit Radius (km)",
        "x_unit": "km",
        "curves": [
            {
                "label": "Total Δv (km/s)",
                "color": "#ff8c3a",
                "fn": lambda r2_km: OrbitalMechanics.hohmann_delta_v(
                    km_to_m(6571), km_to_m(r2_km))["total_dv_ms"] / 1000
            },
            {
                "label": "Δv₁ (km/s)",
                "color": "#3a7fff",
                "fn": lambda r2_km: OrbitalMechanics.hohmann_delta_v(
                    km_to_m(6571), km_to_m(r2_km))["dv1_ms"] / 1000
            },
            {
                "label": "Δv₂ (km/s)",
                "color": "#b06bff",
                "fn": lambda r2_km: OrbitalMechanics.hohmann_delta_v(
                    km_to_m(6571), km_to_m(r2_km))["dv2_ms"] / 1000
            },
        ],
        "x_range": [6800, 42164],
        "annotations": [
            {"x": 42164, "label": "GEO"},
        ]
    },

    # ── PROPULSION ───────────────────────────────────────────────────────────
    "tsiolkovsky": {
        "title": "Tsiolkovsky: Δv vs Mass Ratio",
        "x_label": "Mass Ratio (m₀/mf)",
        "x_unit": "",
        "curves": [
            {
                "label": "Isp 250s (solid)",
                "color": "#ff4a6e",
                "fn": lambda mr: RocketPropulsion.exhaust_velocity(250) * math.log(mr) / 1000
            },
            {
                "label": "Isp 350s (kerosene)",
                "color": "#ff8c3a",
                "fn": lambda mr: RocketPropulsion.exhaust_velocity(350) * math.log(mr) / 1000
            },
            {
                "label": "Isp 450s (hydrolox)",
                "color": "#3a7fff",
                "fn": lambda mr: RocketPropulsion.exhaust_velocity(450) * math.log(mr) / 1000
            },
            {
                "label": "Isp 3000s (ion)",
                "color": "#2affa0",
                "fn": lambda mr: RocketPropulsion.exhaust_velocity(3000) * math.log(mr) / 1000
            },
        ],
        "x_range": [1.1, 20],
        "annotations": [
            {"x": 9.4,  "label": "LEO (~9.4 km/s)"},
        ]
    },

    "twr": {
        "title": "Thrust-to-Weight Ratio vs Thrust",
        "x_label": "Thrust (kN)",
        "x_unit": "kN",
        "curves": [
            {
                "label": "TWR (100t rocket)",
                "color": "#ff8c3a",
                "fn": lambda thrust_kn: RocketPropulsion.twr(thrust_kn * 1000, 100000)
            },
            {
                "label": "TWR (500t rocket)",
                "color": "#3a7fff",
                "fn": lambda thrust_kn: RocketPropulsion.twr(thrust_kn * 1000, 500000)
            },
        ],
        "x_range": [500, 10000],
        "annotations": [
            {"x": 1, "label": "TWR=1 (liftoff)"},
        ]
    },

    # ── PLANETARY ────────────────────────────────────────────────────────────
    "equilibrium_temp": {
        "title": "Planetary Equilibrium Temperature vs Distance",
        "x_label": "Distance from Star (AU)",
        "x_unit": "AU",
        "curves": [
            {
                "label": "Albedo 0.1 (dark)",
                "color": "#ff4a6e",
                "fn": lambda au: PlanetaryScience.equilibrium_temperature(0.1, au)
            },
            {
                "label": "Albedo 0.3 (Earth-like)",
                "color": "#3a7fff",
                "fn": lambda au: PlanetaryScience.equilibrium_temperature(0.3, au)
            },
            {
                "label": "Albedo 0.7 (icy)",
                "color": "#00cfff",
                "fn": lambda au: PlanetaryScience.equilibrium_temperature(0.7, au)
            },
        ],
        "x_range": [0.3, 5.0],
        "annotations": [
            {"x": 0.72, "label": "Venus"},
            {"x": 1.0,  "label": "Earth"},
            {"x": 1.52, "label": "Mars"},
        ],
        "y_ref_lines": [
            {"y": 373, "label": "Boiling point (373K)"},
            {"y": 273, "label": "Freezing point (273K)"},
        ]
    },

    "kepler_period": {
        "title": "Orbital Period vs Semi-major Axis (Kepler III)",
        "x_label": "Semi-major Axis (AU)",
        "x_unit": "AU",
        "curves": [
            {
                "label": "Period (years)",
                "color": "#2affa0",
                "fn": lambda au: seconds_to_days(
                    PlanetaryScience.kepler_third_law_period(au_to_m(au))) / 365.25
            }
        ],
        "x_range": [0.3, 30],
        "annotations": [
            {"x": 0.72,  "label": "Venus"},
            {"x": 1.0,   "label": "Earth"},
            {"x": 1.52,  "label": "Mars"},
            {"x": 5.2,   "label": "Jupiter"},
            {"x": 9.58,  "label": "Saturn"},
            {"x": 19.2,  "label": "Uranus"},
            {"x": 30.05, "label": "Neptune"},
        ]
    },

    # ── RELATIVITY ───────────────────────────────────────────────────────────
    "lorentz_factor": {
        "title": "Lorentz Factor γ vs Velocity",
        "x_label": "Velocity (fraction of c)",
        "x_unit": "c",
        "curves": [
            {
                "label": "γ (Lorentz Factor)",
                "color": "#b06bff",
                "fn": lambda v_c: Relativity.lorentz_factor(v_c * Constants.C)
            },
            {
                "label": "Time Dilation (t/t₀ = γ)",
                "color": "#ff4a6e",
                "fn": lambda v_c: Relativity.lorentz_factor(v_c * Constants.C)
            },
        ],
        "x_range": [0.01, 0.999],
        "annotations": [
            {"x": 0.5,  "label": "0.5c"},
            {"x": 0.9,  "label": "0.9c"},
            {"x": 0.99, "label": "0.99c"},
        ]
    },

    "time_dilation": {
        "title": "Time Dilation: Elapsed Time vs Velocity",
        "x_label": "Velocity (fraction of c)",
        "x_unit": "c",
        "curves": [
            {
                "label": "Observer time (t) for 1yr proper time",
                "color": "#b06bff",
                "fn": lambda v_c: Relativity.time_dilation(
                    365.25 * 24 * 3600, v_c * Constants.C) / (365.25 * 24 * 3600)
            },
        ],
        "x_range": [0.01, 0.999],
        "annotations": [
            {"x": 0.9,  "label": "0.9c → 2.3× slower"},
            {"x": 0.99, "label": "0.99c → 7× slower"},
        ]
    },

    "schwarzschild_radius": {
        "title": "Schwarzschild Radius vs Black Hole Mass",
        "x_label": "Mass (Solar Masses)",
        "x_unit": "M☉",
        "curves": [
            {
                "label": "Schwarzschild Radius (km)",
                "color": "#b06bff",
                "fn": lambda ms: Relativity.schwarzschild_radius(ms * Constants.M_SUN) / 1000
            }
        ],
        "x_range": [1, 100],
        "annotations": [
            {"x": 1,    "label": "1 M☉"},
            {"x": 10,   "label": "Stellar BH (~10 M☉)"},
            {"x": 30,   "label": "GW150914 (~30 M☉)"},
        ]
    },

    "hawking_temperature": {
        "title": "Hawking Temperature vs Black Hole Mass",
        "x_label": "Mass (kg, log scale)",
        "x_unit": "kg",
        "curves": [
            {
                "label": "Hawking Temperature (K)",
                "color": "#ff4a6e",
                "fn": lambda log_m: Relativity.hawking_temperature(10**log_m)
            }
        ],
        "x_range": [10, 35],
        "x_is_log": True,
        "annotations": [
            {"x": math.log10(Constants.M_SUN), "label": "1 M☉"},
        ]
    },

    # ── QUANTUM ───────────────────────────────────────────────────────────────
    "photon_energy": {
        "title": "Photon Energy vs Wavelength",
        "x_label": "Wavelength (nm)",
        "x_unit": "nm",
        "curves": [
            {
                "label": "Energy (eV)",
                "color": "#ffd93a",
                "fn": lambda nm: J_to_eV(
                    QuantumMechanics.photon_energy_from_wavelength(nm * 1e-9))
            }
        ],
        "x_range": [100, 1000],
        "annotations": [
            {"x": 400, "label": "Violet"},
            {"x": 550, "label": "Green"},
            {"x": 700, "label": "Red"},
        ],
        "gradient_x": True
    },

    "wien_law": {
        "title": "Blackbody Spectrum Peak & Radiated Power vs Temperature",
        "x_label": "Temperature (K)",
        "x_unit": "K",
        "curves": [
            {
                "label": "Peak Wavelength (nm)",
                "color": "#ffd93a",
                "fn": lambda T: QuantumMechanics.wien_peak_wavelength(T) * 1e9
            },
            {
                "label": "Power / 1e8 (W/m²)",
                "color": "#ff8c3a",
                "fn": lambda T: QuantumMechanics.stefan_boltzmann(T) / 1e8
            },
        ],
        "x_range": [1000, 30000],
        "annotations": [
            {"x": 3000,  "label": "Red dwarf"},
            {"x": 5778,  "label": "Sun"},
            {"x": 10000, "label": "Blue-white star"},
        ]
    },

    # ── EM ───────────────────────────────────────────────────────────────────
    "coulomb_force": {
        "title": "Coulomb Force vs Separation",
        "x_label": "Separation (mm)",
        "x_unit": "mm",
        "curves": [
            {
                "label": "|F| q=1nC,q=-1nC (N)",
                "color": "#00cfff",
                "fn": lambda mm: abs(ElectroMagnetism.coulomb_force(1e-9, -1e-9, mm * 1e-3))
            },
            {
                "label": "|F| q=10nC,q=-10nC (N)",
                "color": "#3a7fff",
                "fn": lambda mm: abs(ElectroMagnetism.coulomb_force(10e-9, -10e-9, mm * 1e-3))
            },
        ],
        "x_range": [1, 100],
    },
}


@app.route("/api/plot", methods=["POST"])
def get_plot():
    """Generate x/y curve data for a formula plot."""
    data = request.get_json(force=True)
    key  = data.get("formula")

    if key not in PLOT_REGISTRY:
        return jsonify({"ok": False, "error": "No plot available for this formula"}), 404

    p = PLOT_REGISTRY[key]
    x_min, x_max = p["x_range"]
    N = 300

    xs = linspace(x_min, x_max, N)

    curves_out = []
    for curve in p["curves"]:
        ys = []
        for x in xs:
            try:
                y = curve["fn"](x)
                ys.append(y if math.isfinite(y) else None)
            except Exception:
                ys.append(None)
        curves_out.append({
            "label": curve["label"],
            "color": curve["color"],
            "x": xs,
            "y": ys
        })

    return jsonify({
        "ok": True,
        "title": p["title"],
        "x_label": p["x_label"],
        "x_unit": p.get("x_unit", ""),
        "annotations": p.get("annotations", []),
        "y_ref_lines": p.get("y_ref_lines", []),
        "curves": curves_out
    })


@app.route("/api/plots")
def list_plots():
    """Return which formula keys have plots available."""
    return jsonify(list(PLOT_REGISTRY.keys()))


# ─────────────────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    debug = os.environ.get("FLASK_ENV", "development") == "development"
    print(f"\n{'='*50}")
    print("  ASTROCALC DELTA — Mission Control")
    print(f"  http://127.0.0.1:{port}")
    print(f"{'='*50}\n")
    app.run(host="0.0.0.0", port=port, debug=debug)
