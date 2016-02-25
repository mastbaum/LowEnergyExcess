import sys

if len(sys.argv) < 3:
    msg  = '\n'
    msg += "Usage 1: %s output_plot_name $INPUT_ROOT_FILE\n" % sys.argv[0]
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
    if x < 2: continue
    my_proc.add_input_file(sys.argv[x])

# Specify IO mode
my_proc.set_io_mode(fmwk.storage_manager.kREAD)

#Enable filtering
my_proc.enable_filter(True)

# Specify output root file name
my_proc.set_ana_output_file("GenLEERWSourcePlot_out.root")

myfilter = fmwk.MC_LEE_Filter()
my_proc.add_process(myfilter)

plotgen_module = fmwk.GenLEERWSourcePlot()
plotgen_module.setPlotName(sys.argv[1])
my_proc.add_process(plotgen_module)


# Let's run it.
my_proc.run()

# done!
print
print "Finished running ana_processor event loop!"
print

sys.exit(0)

