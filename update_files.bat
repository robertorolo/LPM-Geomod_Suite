@ECHO OFF 
ECHO Updating files from git folder to ar2gems folder
ECHO copying helpers
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\helpers.py" "C:\ar2gems-full-version-python3\Release\Lib" /y

ECHO copying py files
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\py\indicators_transform.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\py\signed_distances_transform.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\py\deterministic.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y

ECHO copying ui files
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\ui\deterministic.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\ui\indicators_transform.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\ui\deterministic.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y


