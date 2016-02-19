import sys, os

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s \'mc\'/\'reco\' $INPUT_ROOT_FILEs $OUTPUT_PATH\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

if sys.argv[1] not in ['reco','mc']:
	msg = '\n'
	msg += 'Specify if you want to use "reco" or "mc" quantities in your first argument to this script!'
	msg += '\n'
	sys.stderr.write(msg)
	sys.exit(1)

from ROOT import gSystem
from ROOT import larlite as fmwk
from ROOT import ertool
from singleE_config import GetERSelectionInstance

# Create ana_processor instance
my_proc = fmwk.ana_processor()
my_proc.enable_filter(True)

use_reco = True if sys.argv[1] == 'reco' else False

# Set input root file
for x in xrange(len(sys.argv)-3):
    my_proc.add_input_file(sys.argv[x+2])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

# Specify output root file name
outfilebase = sys.argv[-1]+'/'+sys.argv[0][:-3]+'_%s'%('mc' if not use_reco else 'reco')
outfile = outfilebase+'.root'
print "%s output file = %s"%(sys.argv[0],outfile)
my_proc.set_ana_output_file(outfile)
my_proc.set_output_file(outfilebase+'_larlite_out.root')

#nueCC beam
eventfilter = fmwk.MC_CCnue_Filter()

LEEana = ertool.ERAnaLowEnergyExcess()
LEEana.SetTreeName("beamNuE")
#LEEana.SetDebug(False)


anaunit = GetERSelectionInstance()
anaunit._mgr.ClearCfgFile()
if not use_reco:
	anaunit._mgr.AddCfgFile(os.environ['LARLITE_USERDEVDIR']+'/SelectionTool/ERTool/dat/ertool_default.cfg')
else:
	anaunit._mgr.AddCfgFile(os.environ['LARLITE_USERDEVDIR']+'/SelectionTool/ERTool/dat/ertool_default_emulated.cfg')


if use_reco:
	anaunit.SetShowerProducer(False,'recoemu')
	anaunit.SetTrackProducer(False,'recoemu')

anaunit._mgr.AddAna(LEEana)
# Add MC filter and analysis unit
# to the process to be run

my_proc.add_process(eventfilter)
my_proc.add_process(anaunit)

#my_proc.run()
my_proc.run()

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)

