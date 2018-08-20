from pathlib2 import Path
import os.path
import socket

this_machine = socket.gethostname()
print("This machine is: " + this_machine)

srcdir = Path(os.path.abspath('../wave_ra/pnnl/'))

tmpdir = Path(os.path.abspath('../wave_ra/tmp/'))

if (this_machine.startswith('lkilcher-') and
        this_machine.endswith('nrel.gov')):
    srcdir = Path('/Volumes/lkilcher/wave_ra/pnnl/')
    tmpdir = Path('/Volumes/lkilcher/wave_ra/tmp/')
