@ECHO OFF 
ECHO Updating files from git folder to ar2gems folder
ECHO copying helpers
xcopy "helpers.py" "C:\ar2gems-full-version-python3\Release\Lib" /y

ECHO copying py files
xcopy "\py\point_transform.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y
xcopy "\py\deterministic.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y
xcopy "\py\vtk_export.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y
xcopy "\py\auto_grid.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y
xcopy "\py\validation.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y
xcopy "\py\block_transform.py" "C:\ar2gems-full-version-python3\Release\plugins\Geostat\python" /y

ECHO copying ui files
xcopy "\ui\point_transform.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y
xcopy "\ui\deterministic.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y
xcopy "\ui\vtk_export.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y
xcopy "\ui\auto_grid.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y
xcopy "\ui\validation.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y
xcopy "\ui\block_transform.ui" "C:\ar2gems-full-version-python3\Release\plugins\Geostat" /y
