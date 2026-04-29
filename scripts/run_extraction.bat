@echo off
echo Running GEE biomass extraction for California...
echo.
"C:\Users\shaht\anaconda3\python.exe" "%~dp0python\01_extract_biomass_gee.py"
echo.
echo Done. Check output\biomass_timeseries.csv
pause
