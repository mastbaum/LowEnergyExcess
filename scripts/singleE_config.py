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

	# Get Default CCSingleE Algorithm instance
	# this information is loaded from:
	# $LARLITE_BASEDIR/python/seltool/GetCCSingleEInstance
	ccsinglee_algo = GetCCSingleEInstance()
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
	Ecut = 50 # in MeV

	anaunit = fmwk.ExampleERSelection()
	anaunit.setDisableXShift(True)

	anaunit._mgr.AddAlgo(pi0_algo)
	anaunit._mgr.AddAlgo(cos_algo)
	anaunit._mgr.AddAlgo(cosmicprimary_algo)
	anaunit._mgr.AddAlgo(cosmicsecondary_algo)
	anaunit._mgr.AddAlgo(cosmicorphanalgo)
	anaunit._mgr.AddAlgo(primary_algo)
	anaunit._mgr.AddAlgo(pid_algo)
	anaunit._mgr.AddAlgo(ccsinglee_algo)
	anaunit._mgr._profile_mode = True

	anaunit.SetMinEDep(Ecut)
	anaunit._mgr._mc_for_ana = True

	return anaunit
