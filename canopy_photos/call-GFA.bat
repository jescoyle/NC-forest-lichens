@echo off
if %1.==.goto usage
if exist %1 (
echo %1 > temp
goto suite )
if exist %1.bmp echo %1.bmp > temp
if not exist %1 if not exist %1.bmp goto filenotfound
:suite
echo 1 >> temp
echo 936 0 >> temp
echo 0 936 >> temp
echo 1872 936 >> temp
echo %2 >> temp
echo n >> temp
echo 18 24 >> temp
echo %3 >> temp
GFA-WIN < temp
del gf_%1.txt
ren gapfract.txt gf_%1.txt
del temp
goto end

:usage
echo USAGE : call-GFA <file_name> <c> <site_declination>
goto end

:filenotfound
echo file not found
goto end

:end
