import subprocess
import sys

pkgs = ['scipy', 'pyevtk', 'matplotlib', 'seaborn', 'sklearn']
for pkg in pkgs:
    subprocess.call([sys.executable, '-m', 'pip', 'install', '{}'.format(pkg)])
