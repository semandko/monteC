@echo off
rem Set the working directory to the location of your C++ program
cd c:\git_rep\GitHub\isp\monteC\

rem Profile your program
gprof static.exe gmon.out > analysis.txt

rem Pause the script to keep the command prompt window open (optional)
pause
