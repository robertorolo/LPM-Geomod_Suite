import subprocess
import sys

pkgs = ['scipy', 'pyevtk', 'matplotlib']
for pkg in pkgs:
    subprocess.call([sys.executable, '-m', 'pip', 'install', '{}'.format(pkg)])
