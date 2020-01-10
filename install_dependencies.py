import subprocess
import sys

pkg_lst = ['scipy', 'pyevtk']
for pkg in pkg_lst:
    subprocess.call([sys.executable, '-m', 'pip', 'install', '{}'.format(pkg)])
