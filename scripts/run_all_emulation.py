import os, datetime, sys

empath = '/Users/davidkaleko/larlite/UserDev/LArLiteApp/RecoEmulatorApp/mac/'
em = 'run_reco_emulator.py'
cfg = 'recoemu_master.fcl'
outbase = '/Users/davidkaleko/Data/larlite/joseph_LEE_files/'
smps = [ 'osc_bnb_70kv_all_', 'osc_cosmics_70kv_all_' ]

for sample in smps:
	cmd = 'python %s %s %s %s'%(empath+em,empath+cfg,outbase+sample+'mcinfo.root',outbase+sample+'recoemmaster.root')
	print "emulation command: %s"%cmd
	os.system(cmd)
