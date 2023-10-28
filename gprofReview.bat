@echo off
rem Set the working directory to the location of your C++ program
cd C:\devc++

rem Profile your program
gprof Project1.exe gmon.out > analysis.txt

rem Pause the script to keep the command prompt window open (optional)
pause
