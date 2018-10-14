from pathlib2 import Path
import os.path
import socket

this_machine = socket.gethostname()
print("This machine is: " + this_machine)

srcdir = Path(os.path.abspath('../wave_ra/pnnl/'))

tmpdir = Path(os.path.abspath('./tmpdata/'))

projdir = Path(os.path.abspath(__file__)).parent

if this_machine.startswith('lkilcher-26339s') or \
   this_machine.startswith('lkilcher-32045s'):
    srcdir = Path(os.path.expanduser('~/tmp/wave_ra/pnnl/'))

# elif (this_machine.startswith('lkilcher-32045s') and
#         this_machine.endswith('nrel.gov')):
#     srcdir = Path('/Volumes/lkilcher/wave_ra/pnnl/')
