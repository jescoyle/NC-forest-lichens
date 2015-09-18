if exist %1 goto runprogram
if not exist %1 goto filenotfound
:runprogram
SET "DECL="
SET "LON="
SET "LAT="
SET "ELEV="
FOR /F "tokens=1,2,3,4 delims=," %%A IN (%1) DO (
		echo Parameter file is %1
		SET DECL=%%A
		SET LON=%%B
		SET LAT=%%C
		SET ELEV=%%D
)
(
		echo 0
		echo 0
		echo 0
		echo %ELEV%
		echo %LAT%
		echo %LON%
		echo 5
		echo 10
		echo 0.1
		echo 0.11
		echo 0.16
		echo 0.033
		echo 0.18
		echo 12
		echo 1 32 60 91 121 152 182 213 244 274 305 335
		echo 0.0
		echo 0
		echo 0
	) > parms.txt
FOR %%G IN (*.bmp) DO (
	
	call-GFA.bat %%~nG c %DECL%
	
	PARCLR parms.txt gf_%%~nG.txt rad_%%~nG.txt
)
goto end

:filenotfound
echo parameter file not found

:end
