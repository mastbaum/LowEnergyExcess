import sys, os

if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILEs\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from ROOT import gSystem
from ROOT import larlite as fmwk
from ROOT import ertool

# Create ana_processor instance
my_proc = fmwk.ana_processor()
my_proc.enable_filter(True)

# Set input root file
for x in xrange(len(sys.argv)):
    if not x: continue
    my_proc.add_input_file(sys.argv[x])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

my_proc.set_ana_output_file('crycorsikadebug.root')

mod_debug = ertool.ERAnaCryCorsikaDebug()

anaunit = fmwk.ExampleERSelection()
anaunit.SetShowerProducer(True,'mcreco')
anaunit.SetTrackProducer(True,'mcreco')
anaunit.SetMinEDep(10.)
anaunit.setDisableXShift(True)
anaunit._mgr._mc_for_ana = True

anaunit._mgr.ClearCfgFile()
anaunit._mgr.AddCfgFile(os.environ['LARLITE_USERDEVDIR']+'/SelectionTool/ERTool/dat/ertool_default.cfg')
anaunit._mgr.AddAna(mod_debug)
# Add MC filter and analysis unit
# to the process to be run
my_proc.add_process(anaunit)

my_proc.run(0,1000)


# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)

