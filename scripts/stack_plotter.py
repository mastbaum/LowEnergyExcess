#####################################################################
# Please do not modify this and push it to the repository unless all other uses of this analysis
# are aware of what you are changing. If you want to make changes to this script, please do it 
# in your own local area (IE in an ipython notebook or something).
#####################################################################

# This is the production script to make stacked backgrounds for the Low Energy Excess analysis.
# This script assumes you have already run the various "singleE_*_selection.py" scripts to actually
# run the analysis and generate ttrees.
# All this script does is read in those TTrees, apply analysis cuts, and make the stacked
# background histogram from the TTrees

# To use this script, you will need to have the following packages installed:
# numpy
# matplotlib
# pandas

# To figure out how to install these things, google it. If you use pip, consider
# >> pip install numpy
# for example.

# Also in this script are hard-coded paths to where the input files live.
# If yours live in a different place, you will have to modify this.

#from ROOT import TFile, TTree
import os
import numpy as np
import matplotlib.pylab as plt
import pandas as pd
from root_numpy import root2array
from collections import OrderedDict

#####################################################################
##  These are some sample analysis cuts one may want to apply
##  Look at how these are formatted and used if you want to apply your own
#10cm from all sides
defaultcut = '_e_Edep > 60.'
fidxmin, fidxmax = 10., 246.
fidymin, fidymax = -106.5, 86.5
fidzmin, fidzmax = 10., 1026.8
fidvolcut = '_y_vtx > %f and _y_vtx < %f '%(fidymin, fidymax)
fidvolcut += 'and _z_vtx > %f and _z_vtx < %f'%(fidzmin, fidzmax)
fidvolcut += ' and _x_vtx > %f and _x_vtx < %f'%(fidxmin, fidxmax)
BGWstart = 3.5
BGWend   = 5.2
BGWcut = '_flash_time > %f and _flash_time < %f'%(BGWstart,BGWend)
#####################################################################

# Default x-axis variable for stacked histograms
default_plot_variable = '_e_nuReco_better'

# Where the output files live that contain ttrees to plot from
#filebase = os.environ['LARLITE_USERDEVDIR']+'/LowEnergyExcess/output/'
filebase = '/Users/davidkaleko/Data/larlite/nevis_LEE_results/'

# Whether you use perfect reco input or reco emulated input
mc_or_reco = 'mc'

filenames = OrderedDict([
('nue','nue/%s_alljobs_TNfinal.root'%(mc_or_reco)),
('numu','numu/%s_alljobs_TNfinal.root'%(mc_or_reco)),
('nc','nc/%s_alljobs_TNfinal.root'%(mc_or_reco)),
('bite','dirt/%s_alljobs_TNfinal.root'%(mc_or_reco)),
('cosmic','cosmic/%s_alljobs_TNfinal.root'%(mc_or_reco)),
('cosmicoutoftime','cosmicoutoftime/%s_alljobs_TNfinal.root'%(mc_or_reco)),
('lee','lee/%s_alljobs_TNfinal.root'%(mc_or_reco))])

treenames = { 'nue' : 'beamNuE',
             'cosmic' : 'cosmicShowers',
             'cosmicoutoftime':'cosmicOutOfTime',
             'numu' : 'beamNuMu',
             'nc' : 'beamNC',
             'bite': 'dirt',
            'lee' : 'LEETree'}

labels = { 'nue' : 'Beam Intrinsic $\\nu_e$',
          'cosmicoutoftime':'Cosmic Out Of Time, Nu In Time',
         'cosmic' : 'Corsika Cosmics In Time (BGW Exposure Scaled)',#'CRY Cosmic, in-time',
         'numu' : 'Beam Intrinsic $\\nu_\mu$',
         'nc' : 'Beam Intrinsic NC', 
          'bite': 'B.I.T.E. (In Cryostat)',
         'lee' : 'Scaled Signal'}

colors = { 'nue' : '#269729', #kGreen-2
         'numu' : '#4B4EAC', #kBlue-5
          'nc' : '#6B70F5', #kBlue-9
          'cosmic' : '#D12C2C', #kRed-3
          'cosmicoutoftime' : '#F700FF', #temporary pink
          'bite' : '#9F9F9F', #grayish?
          'lee' : '#E65C00' #orangish
          }


# These weights are to scale to 6.6e20 POT
# You may need to change the number of events generated
# (how many events you ran over BEFORE any event filtering)
# Keep the 'lee' weight at 1, but in the singleE_lee_selection
# script, there you will need to input the number of events
# you run over (AFTER the event filtering).
# You may need to run over the lee sample once just to see
# what number to put into that run script.
scaling_weights = \
{ 'nue' : 6.6e20/(3.1845e17*35880),
 #(211,000 ms total exposure)/(7.2ms * 36600 evts generated)
  'cosmic' : 211000/(7.25*36600), 
   #'cosmic' : (1.056e6)/(36600),
  'numu' : 6.6e20/((1.203e15)*98860),
  'nc' : 6.6e20/(1.203e15*98860),
  'bite' : 6.6e20/(1.203e15*98860),
  'cosmicoutoftime': 6.6e20/(1.202e15*49820),
  'lee' : 1 }


# Read in all the ttrees to pandas dataframes
dfs = OrderedDict()
for key, filename in filenames.iteritems():
    dfs.update( { key : pd.DataFrame( root2array( filebase + filename, treenames[key] ) ) } )

#throw away intime cosmics from outoftime sample
dfs['cosmicoutoftime']=dfs['cosmicoutoftime'].query('_mc_time<3100 or _mc_time>4700')
#enforce the out of time cosmics truly come from cosmics
dfs['cosmicoutoftime']=dfs['cosmicoutoftime'].query('_mc_origin == 2')

#Hack the open cosmics so they all flash in the middle of the BGW
#(this is how we scale cosmics to total BGW exposure time..
# this way we can just apply the same flashmatch cut and not have
# to do any other gymnastics specific to this one sample)
dfs['cosmic']['_flash_time'] = ((BGWstart+BGWend)/2.)

# Uncomment this if you want to see what variables are stored in the dataframes
#dfs['cosmic'].info()

# This function makes a weighted numpy histogram from the dataframes
# You give it a "query" (analysis cut), the variable you want to have
# on the x-axis, and an optional "scale factor" which is just to convert
# units (IE from MEV to GEV)
# It returns a dictionary of { sample name : histogram }
def gen_histos( binning = np.linspace(0,10,1), myquery='', plotvar = default_plot_variable, scalefactor = 1.):
    nphistos = OrderedDict()

    for key, df in dfs.iteritems():
      mydf = df.query(myquery) if myquery else df
      # Cosmic weights are applied slightly differently
      # Once we switch to full flash-matching (and not just
      # scaling to beam gate open exposure time) this will
	    # change.
      if key == 'cosmic':
          myweights = np.ones(mydf[default_plot_variable].shape[0])
      else:
          myweights = np.array(mydf['_weight'])
      myweights *= scaling_weights[key]
      nphistos.update( {key : np.histogram(mydf[plotvar]/scalefactor,
                                   bins=binning,
                                   weights=myweights)} )
    return nphistos

# This function loops over the dataframes and makes the 
# stacked background. It takes as input the dictionary of already-created
# histograms and just draws them prettily
def plot_fullstack( binning = np.linspace(0,10,1), myquery='', plotvar = default_plot_variable, \
                    scalefactor = 1., user_ylim = None):

    fig = plt.figure(figsize=(10,6))
    plt.grid(True)
    lasthist = 0
    myhistos = gen_histos(binning=binning,myquery=myquery,plotvar=plotvar,scalefactor=scalefactor)
    for key, (hist, bins) in myhistos.iteritems():

      plt.bar(bins[:-1],hist,
              width=bins[1]-bins[0],
              color=colors[key],
              bottom = lasthist,
              edgecolor = 'k',
              label='%s: %d Events'%(labels[key],sum(hist)))
      lasthist += hist
     

    plt.title('CCSingleE Stacked Backgrounds',fontsize=25)
    plt.ylabel('Events',fontsize=20)
    if plotvar == '_e_nuReco' or plotvar == '_e_nuReco_better':
        xstring = 'Reconstructed Neutrino Energy [GeV]' 
    elif plotvar == '_e_CCQE':
        xstring = 'CCQE Energy [GeV]'
    else:
        xstring = plotvar
    plt.xlabel(xstring,fontsize=20)
    plt.legend()
    plt.xticks(list(plt.xticks()[0]) + [binning[0]])
    plt.xlim([binning[0],binning[-1]])
    #plt.ylim([0,130])

if __name__ == '__main__':
  mybins = np.linspace(0.1,3.0,15)
  mycuts = defaultcut + ' and ' + BGWcut + ' and ' + fidvolcut
  print "The cuts used to make the plot you're seeing are: "
  print mycuts
  plot_fullstack(binning=mybins, myquery=mycuts, \
    plotvar='_e_nuReco', scalefactor=1000.)
  plt.show()


