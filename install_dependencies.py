import subprocess
import sys
packs = ['scipy', 'pyevtk']
for pkg in packs:
    subprocess.call([sys.executable, '-m', 'pip', 'install', '{}'.format(pkg)])