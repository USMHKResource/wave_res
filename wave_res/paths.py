from pathlib2 import Path
import os
import socket

this_machine = socket.gethostname()
print("This machine is: " + this_machine)

srcdir = Path(os.path.abspath('../wave_ra/pnnl/'))

tmpdir = Path(os.path.abspath('./tmpdata/'))

pkgdir = Path(os.path.abspath(__file__)).parent

if this_machine.startswith('lkilcher-26339s') or \
   this_machine.startswith('lkilcher-32045s'):
    srcdir = Path('/Volumes/lkilcher/wave_ra/pnnl/')
    maskdir = Path(os.path.expanduser('~/Dropbox/tmp/wave_ra/icemask/'))
    #srcdir = Path(os.path.expanduser('~/tmp/wave_ra/pnnl/'))

# Gabriel's Paths
if this_machine.startswith('constance') or this_machine.startswith('node'):
    rootFld = Path('/pic/projects/fvwetland/')
    srcdir = rootFld / 'gabriel/waveEnergyResource/assessment/hindcast/'
    tmpdir = srcdir / 'resource/tmpFreq/'

#Aidans Paths
if this_machine.startswith('abharath-34229s'):
    srcdir = Path('/frequencyResults/')
    # 'Y:\wind\public\users\lkilcher\wave_ra\pnnl\'

if 'maskdir' not in vars():
    maskdir = srcdir / 'iceMask/' / 'akMask'


def mkdir(directory):
    if not os.path.exists(directory):
        os.makedirs(directory)
    else:
        print('directory exists: ' + directory)
