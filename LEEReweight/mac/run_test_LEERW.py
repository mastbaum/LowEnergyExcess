import sys



if len(sys.argv) < 2:
    msg  = '\n'
    msg += "Usage 1: %s $INPUT_ROOT_FILE\n" % sys.argv[0]
    msg += '\n'
    sys.stderr.write(msg)
    sys.exit(1)

from ROOT import gSystem
from ROOT import larlite as fmwk

# Create ana_processor instance
my_proc = fmwk.ana_processor()
#my_proc.set_verbosity(fmwk.msg.kNORMAL)

# Set input root file
for x in xrange(len(sys.argv)):
    if x == 0: continue
    my_proc.add_input_file(sys.argv[x])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

# Specify output root file name
my_proc.set_ana_output_file("test_LEERW_anaout.root")
my_proc.enable_filter(True)
my_proc.add_process(fmwk.MC_LEE_Filter())
test_module = fmwk.test_LEERW()
test_module.set_debug(False)
my_proc.add_process(test_module)

# Let's run it.
my_proc.run()#1,50

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)

