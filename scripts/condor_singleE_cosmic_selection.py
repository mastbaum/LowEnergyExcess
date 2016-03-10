import sys, os

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: input_config_file.txt\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from ROOT import gSystem
from ROOT import larlite as fmwk
from ROOT import ertool
from singleE_config import GetERSelectionInstance

# Create ana_processor instance
my_proc = fmwk.ana_processor()
#my_proc.enable_event_alignment(False)
my_proc.enable_filter(True)

emulator_cfg = ''

#Parse the config file
with open(sys.argv[1]) as f:
    content = f.readlines()
newcontent = [x.strip('\n') for x in content]
use_reco = True if newcontent[0] == 'reco' else False
if use_reco:
    emulator_cfg = newcontent[1]
    newcontent.pop(1)

# Set input root files
for infile in newcontent[1:-1]:
    my_proc.add_input_file(infile)

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kBOTH)

# Specify output root file name
outfile = newcontent[-1]
print "%s output file = %s"%(sys.argv[0],outfile)
my_proc.set_ana_output_file(outfile)
my_proc.set_output_file(outfile[:-5]+'_larlite_out.root')

#nueCC beam
eventfilter = fmwk.MC_cosmic_Filter()

LEEana = ertool.ERAnaLowEnergyExcess()
LEEana.SetTreeName("cosmicShowers")
#LEEana.SetDebug(False)

anaunit = GetERSelectionInstance()
anaunit._mgr.ClearCfgFile()

#DISBABLE x-shift if you are using open cosmics and scaling to total BGW exposure
anaunit.setDisableXShift(True)

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
#Add reco emulator if necessary!
if use_reco:
    emulator = fmwk.EmuDriver()
    emulator.set_config(emulator_cfg)
    print "USING RECO EMULATOR. CONFIG FILE USED FOR EMULATOR IS %s"%emulator_cfg
    my_proc.add_process(emulator)

my_proc.add_process(anaunit)

my_proc.run()

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)

