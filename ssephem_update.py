import os
import sys
baryc_dir = os.path.dirname(os.path.abspath(__file__))+'/SSEphem/'
#baryc_dir = '/home/ronny/SSEphem/'          # Testing
sys.path.append(baryc_dir)
import update_ssephem

os.chdir(baryc_dir)

update_ssephem.SSEphemDownload()           # DEc403
#update_ssephem.SSEphemDownload_405()        # DEc405
update_ssephem.LeapSecUpdate()
update_ssephem.IersUpdate()
