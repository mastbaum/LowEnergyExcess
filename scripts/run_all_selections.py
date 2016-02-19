import os, datetime, sys

_use_reco = False

input_base = '/Users/davidkaleko/Data/larlite/joseph_LEE_files/'
output_dir = '/Users/davidkaleko/larlite/UserDev/LowEnergyExcess/output/70KV/crorphanfix/mcinfo_only/'#mc_startdir/'
if _use_reco:
	output_dir = '/Users/davidkaleko/larlite/UserDev/LowEnergyExcess/output/70KV/recoemmaster/with_reco_noenergysmear_retrained/'


cosmics_files =  input_base + 'osc_cosmics_70kv_all_mcinfo.root'
cosmics_files += ' ' + input_base + 'osc_cosmics_70kv_all_opdata.root'

if _use_reco:
	# cosmics_files += ' ' + input_base + 'osc_cosmics_70kv_all_reco3d_fuzzyshower.root'
	# cosmics_files += ' ' + input_base + 'osc_cosmics_70kv_all_reco3d_kalmanhitcc.root'
	cosmics_files += ' ' + input_base + 'osc_cosmics_70kv_all_recoem_master_noshowereneegysmear.root'

dirt_files = '/Users/davidkaleko/Data/larlite/mcc7/mcc7_test_bnb_nu_v1/larlite_mcc7_test_bnb_nu_v1_mcinfo_all.root'
dirt_files += ' /Users/davidkaleko/Data/larlite/mcc7/mcc7_test_bnb_nu_v1/larlite_mcc7_test_bnb_nu_v1_opreco_all.root'

bnb_files =  input_base + 'osc_bnb_70kv_all_mcinfo.root'
bnb_files += ' ' + input_base + 'osc_bnb_70kv_all_opdata.root'

if _use_reco:
	# bnb_files += ' ' + input_base + 'osc_bnb_70kv_all_reco3d_fuzzyshower.root'
	# bnb_files += ' ' + input_base + 'osc_bnb_70kv_all_reco3d_kalmanhitcc.root'
	bnb_files += ' ' + input_base + 'osc_bnb_70kv_all_recoem_master_noshowereneegysmear.root'
		
lee_files = '/Users/davidkaleko/Data/larlite/LEE_generation/LEEgen_mcinfo_all.root'

if _use_reco:
	lee_files += ' /Users/davidkaleko/Data/larlite/LEE_generation/LEEgen_mcinfo_recoemtest.root'

starttime = datetime.datetime.now()
print "run_all_selections start time is",starttime
if cosmics_files: os.system('python singleE_cosmic_selection.py %s %s %s'%('reco' if _use_reco else 'mc',cosmics_files,output_dir))
if dirt_files: os.system('python singleE_dirt_selection.py %s %s %s'%('reco' if _use_reco else 'mc',dirt_files,output_dir))
if bnb_files: os.system('python singleE_nc_selection.py %s %s %s'%('reco' if _use_reco else 'mc',bnb_files,output_dir))
if bnb_files: os.system('python singleE_nue_selection.py %s %s %s'%('reco' if _use_reco else 'mc',bnb_files,output_dir))
if bnb_files: os.system('python singleE_numu_selection.py %s %s %s'%('reco' if _use_reco else 'mc',bnb_files,output_dir))
# if lee_files: os.system('python singleE_LEE_selection.py %s %s %s'%('reco' if _use_reco else 'mc',lee_files,output_dir))
print "run_all_selections total time duration is",datetime.datetime.now()-starttime
