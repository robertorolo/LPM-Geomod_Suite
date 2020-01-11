import subprocess
import sys

pkg = 'scipy pyevtk'
subprocess.call([sys.executable, '-m', 'pip', 'install', '{}'.format(pkg)])
