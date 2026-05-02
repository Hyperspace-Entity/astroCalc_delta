@echo off
cd /d "%~dp0"

echo.
echo  =============================================
echo   AstroCalc Delta -- Local Setup
echo  =============================================
echo.

python --version >nul 2>&1
IF ERRORLEVEL 1 (
    echo  [ERROR] Python not found. Install from https://python.org
    pause
    exit /b 1
)

echo  [1/3] Creating virtual environment...
python -m venv venv

echo  [2/3] Installing dependencies...
call venv\Scripts\activate.bat
pip install flask flask-cors gunicorn --quiet
pip install flask-cors --quiet

echo  [3/3] Launching AstroCalc Delta...
echo.
echo  Open your browser at: http://127.0.0.1:5000
echo  Press Ctrl+C to stop the server.
echo.

python app.py

pause
