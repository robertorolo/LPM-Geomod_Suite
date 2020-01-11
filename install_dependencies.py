import subprocess
import sys

pkg = 'scipy'
subprocess.call([sys.executable, '-m', 'pip', 'install', '{}'.format(pkg)])
