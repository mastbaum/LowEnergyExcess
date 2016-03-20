from ROOT import gSystem
from ROOT import larlite as fmwk
from ROOT import ertool
from seltool.ccsingleeDef import GetCCSingleEInstance
from seltool.primaryfinderDef import GetPrimaryFinderInstance
from seltool.trackpidDef import GetTrackPidInstance
#from seltool.cosmictaggerDef import GetCosmicTaggerInstance
from seltool.trackDresserDef import GetTrackDresserInstance
from seltool.primarycosmicDef import GetPrimaryCosmicFinderInstance
from seltool.pi0algDef import GetERAlgoPi0Instance

def GetERSelectionInstance():

	# Make an instance of ERAlgoFlashMatch using defaults defined in ertool_default(_mc).cfg
	flashmatch_algo = ertool.ERAlgoFlashMatch()
	flashmatch_algo.SetIgnoreShowers(False)
	flashmatch_algo.SetIgnoreCosmics(True)

	# Get Default CCSingleE Algorithm instance
	# this information is loaded from:
	# $LARLITE_BASEDIR/python/seltool/GetCCSingleEInstance
	ccsinglee_algo = GetCCSingleEInstance()
	#ccsinglee_algo.SetVerbosity(0)
	#ccsinglee_algo.setVerbose(False)

	# primary finder algorithm
	# this information is loaded from:
	# $LARLITE_BASEDIR/python/seltool/GetPrimaryFinderInstance
	primary_algo = GetPrimaryFinderInstance()

	# primary cosmic algoithm 
	# this information is loaded from:
	# $LARLITE_BASEDIR/python/seltool/primarycosmicDef.py
	cosmicprimary_algo = GetPrimaryCosmicFinderInstance()
	cosmicsecondary_algo = ertool.ERAlgoCRSecondary()
	cosmicorphanalgo = ertool.ERAlgoCROrphan()
	# track PID algorithm
	# this information is loaded from:
	# $LARLITE_BASEDIR/python/seltool/GetTrackPidInstance
	pid_algo = GetTrackPidInstance()
	#pid_algo.setVerbose(False)

	# cosmic tagger algo
	#cos_algo = GetCosmicTaggerInstance()
	cos_algo = GetTrackDresserInstance()
	#cos_algo.setVerbose(False)

	pi0_algo = GetERAlgoPi0Instance()
	# here set E-cut for Helper & Ana modules
	#This cut is applied in helper... ertool showers are not made if the energy of mcshower or reco shower
	#is below this threshold. This has to be above 0 or else the code may segfault. This is not a "physics cut".
	#Do not change this value unless you know what you are doing.
	#Ecut = 50 # in MeV
	Ecut = 10 #temporary trying this to see if it helps pi0 mids at low energy
	
	#anaunit = fmwk.ERSelSaveSingleEEvents()
	anaunit = fmwk.ExampleERSelection()
	anaunit.SetShowerProducer(True,'mcreco')
	anaunit.SetTrackProducer(True,'mcreco')

	anaunit.SetFlashProducer('opflashSat')

	anaunit.setDisableXShift(False)

	anaunit._mgr.AddAlgo(ertool.ERAlgoTagEmulatedDeletionsCosmic())
	
	# pi0 algo takes a long time on cosmics files (may showers)...
	# first run track dresser to gobble up most of the showers and 
	# pi0 algo will run much faster (I hope!)
	anaunit._mgr.AddAlgo(cos_algo)
	anaunit._mgr.AddAlgo(pi0_algo)
	
	anaunit._mgr.AddAlgo(cosmicprimary_algo)
	anaunit._mgr.AddAlgo(cosmicsecondary_algo)
	anaunit._mgr.AddAlgo(cosmicorphanalgo)
	anaunit._mgr.AddAlgo(primary_algo)
	# anaunit._mgr.AddAlgo(pid_algo)
	anaunit._mgr.AddAlgo(ccsinglee_algo)
	# Is this where flashmatch_algo should go?
	# First we reconstruct nues and all that, then say if the electron's associated flash
	# is outside of the BGW we throw it out?
	# IE that's an analysis cut later by asking about flash_time?
	# Or maybe it should go right before ccsinglee_algo?
	# That way when it tags something as kcosmic, algosingleE auto-ignores it?
	# However, right now all solo shower particles are tagged as cosmic by the flash-matcher
	# Because it only works for tracks! For now, flashmatch_algo has to live after
	# and it has to be an analysis cut. This will be changed ASAP.
	anaunit._mgr.AddAlgo(flashmatch_algo)

	# Testing adding this... it looks for flashes shared b/t the neutrino and others
	# and potentially adds the "others" as children of the neutrino, or tags
	# the neutrino as a pi0 MID
	#anaunit._mgr.AddAlgo(ertool.ERAlgoNueSharedFlashMerger())
	anaunit._mgr._profile_mode = True

	anaunit.SetMinEDep(Ecut)
	anaunit._mgr._mc_for_ana = True

	return anaunit
