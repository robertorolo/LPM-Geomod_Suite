@ECHO OFF 
ECHO Updating files from git folder to ar2gems folder
ECHO copying helpers
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\helpers.py" "C:\ar2gems-full-version-python3\Release\Lib" /y

ECHO copying py files
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\py\transform.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\py\deterministic.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\py\vtk_export.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\py\auto_grid.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\py\validation.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\py\block_transform.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y

ECHO copying ui files
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\ui\transform.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\ui\deterministic.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\ui\vtk_export.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\ui\auto_grid.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\ui\validation.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y
xcopy "C:\Users\Roberto Rolo\Documents\GitHub\LPM-Geomod_Suite\ui\block_transform.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y