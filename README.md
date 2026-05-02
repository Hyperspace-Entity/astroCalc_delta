# AstroCalc Delta

A production-grade physics calculator covering 6 domains and 30+ equations.
Built with Flask (Python) + vanilla JS. Runs locally, deploys to Railway in minutes.

---

## Domains

| Domain | Formulas |
|---|---|
| 🪐 Orbital Mechanics | Orbital velocity, period, escape velocity, Hohmann transfer, vis-viva, eccentricity, inclination Δv, GEO altitude |
| 🚀 Rocket Propulsion | Tsiolkovsky, mass ratio, propellant fraction, thrust, TWR, burn time, Isp |
| 🌍 Planetary Science | Surface gravity, Kepler III, Roche limit, equilibrium temperature, synodic period, mean density |
| ⚫ Relativity & BH | Lorentz factor, time dilation, length contraction, E=mc², Schwarzschild radius, Hawking temperature, Kerr horizon, gravitational redshift |
| ⚡ Electromagnetism | Coulomb's law, cyclotron radius/frequency, Larmor radiation |
| ⚛️ Quantum Mechanics | Photon energy, de Broglie wavelength, hydrogen levels, Heisenberg uncertainty, Wien's law, Compton scattering |

---

## Local Setup (Windows)

### Option A — Double-click (easiest)
1. Double-click `run.bat`
2. Open `http://127.0.0.1:5000` in your browser

### Option B — Manual
```bat
python -m venv venv
venv\Scripts\activate
pip install -r requirements.txt
python app.py
```

---

## Project Structure

```
astroCalc_delta/
├── app.py            ← Flask server + API routes + formula registry
├── math_module.py    ← Pure physics engine (no dependencies)
├── index.html        ← Single-file frontend
├── requirements.txt  ← Flask, flask-cors, gunicorn
├── Procfile          ← Railway / Render start command
├── railway.json      ← Railway config
├── nixpacks.toml     ← Build config for Railway
├── run.bat           ← Windows one-click start
└── .gitignore
```

---

## API Reference

### `GET /api/formulas`
Returns the full formula catalogue with labels, equations, and parameter definitions.

### `POST /api/calculate`
```json
{
  "formula": "tsiolkovsky",
  "params": {
    "isp": 450,
    "m0_kg": 10000,
    "mf_kg": 3000
  }
}
```
Response:
```json
{
  "ok": true,
  "results": {
    "Delta-v": "6722.45 m/s",
    "Delta-v (km/s)": "6.7224 km/s"
  }
}
```

### `GET /api/constants`
Returns key physical constants with units.

### `GET /api/health`
Health check endpoint used by Railway.

---

## Deploy to Railway

1. Push this folder to a GitHub repository
2. Go to [railway.app](https://railway.app) → New Project → Deploy from GitHub
3. Select your repo
4. Railway auto-detects Python via `nixpacks.toml` and uses the `Procfile`
5. Add environment variable if needed: `FLASK_ENV=production`
6. Done — Railway gives you a public URL

---

## Adding a New Formula

In `app.py`, add an entry to `FORMULA_REGISTRY`:

```python
"my_formula": {
    "label": "My Formula Name",
    "domain": "orbital",   # orbital | propulsion | planetary | relativity | em | quantum
    "equation": "x = a·b",
    "params": [
        {"id": "a", "label": "Value A", "unit": "m", "default": 1.0},
        {"id": "b", "label": "Value B", "unit": "s", "default": 1.0},
    ],
    "fn": lambda p: {
        "Result": f"{p['a'] * p['b']:.4f} m/s"
    }
},
```

The frontend picks it up automatically — no JS changes needed.
