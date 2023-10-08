@echo off
rem Check if ImageMagick is installed
magick -version > nul 2>&1
if %errorlevel% neq 0 (
    echo ImageMagick is not installed. Please install it and make sure 'magick' is in your PATH.
    exit /b 1
)

rem Convert PPM to JPEG
magick ukr.ppm ukr.jpeg

rem Check if conversion was successful
if %errorlevel% neq 0 (
    echo Conversion failed.
    exit /b 1
)

echo Conversion completed successfully.
exit /b 0
