#!/usr/bin/env python

'''Run all event selections.

Based on run_all_selections.py by D. Kaleko.

A. Mastbaum <mastbaum@uchicago.edu>, 2016/11/16
'''

import argparse
import datetime
import os

from singleE_selection import run

parser = argparse.ArgumentParser(description='Run all sample selection')
parser.add_argument('-r', '--reco', action='store_true',
                    help='Use reconstructed quantities')
parser.add_argument('out_path', default='.', help='Output file path')
parser.add_argument('in_path', default='.', help='Input file path')
args = parser.parse_args()

dev_path = os.getenv('LARLITE_USERDEVDIR')
script_base = os.path.join(dev_path, 'LowEnergyExcess', 'scripts')
script = os.path.join(script_base, 'singleE_%s_selection.py')

signals = {
    'cosmic': [
        os.path.join(args.in_path, 'scan_prodcosmics_corsika_cmc_uboone_mcc7_detsim_v2.root'),
        os.path.join(args.in_path, 'scan_prodcosmics_corsika_cmc_uboone_mcc7_detsim_v2_recoemmaster_out.root'),
    ],
    'cosmicoutoftime': [
        os.path.join(args.in_path, 'scan_prodgenie_bnb_nu_cosmic_uboone_mcc7_detsim_v1.root'),
        os.path.join(args.in_path, 'scan_prodgenie_bnb_nu_cosmic_uboone_mcc7_detsim_v1_recoemmaster_out.root'),
    ],
    'nue': [
        os.path.join(args.in_path, 'scan_prodgenie_bnb_intrinsic_nue_uboone_mcc7_detsim_v2.root'),
        os.path.join(args.in_path, 'scan_prodgenie_bnb_intrinsic_nue_uboone_mcc7_detsim_v2_recoemmaster_out.root'),
    ],
    'numu': [
        os.path.join(args.in_path, 'scan_prodgenie_bnb_nu_uboone_mcc7_detsim_v2.root'),
        os.path.join(args.in_path, 'scan_prodgenie_bnb_nu_uboone_mcc7_detsim_v2_recoemmaster_out.root'),
    ],
    'nc': [
        os.path.join(args.in_path, 'scan_prodgenie_bnb_nu_uboone_mcc7_detsim_v2.root'),
        os.path.join(args.in_path, 'scan_prodgenie_bnb_nu_uboone_mcc7_detsim_v2_recoemmaster_out.root'),
    ],
    'dirt': [
        os.path.join(args.in_path, 'scan_prodgenie_bnb_nu_uboone_mcc7_detsim_v2.root'),
        os.path.join(args.in_path, 'scan_prodgenie_bnb_nu_uboone_mcc7_detsim_v2_recoemmaster_out.root'),
    ],
    'LEE': [
        os.path.join(args.in_path, 'scan_prodgenie_bnb_intrinsic_nue_uboone_mcc7_detsim_v2.root'),
        os.path.join(args.in_path, 'scan_prodgenie_bnb_intrinsic_nue_uboone_mcc7_detsim_v2_recoemmaster_out.root'),
    ],
}

starttime = datetime.datetime.now()

for k, v in signals.items():
    print '== Selection:', k, '=='
    if not args.reco:
        v = v[:1]
    run(k, v, args.out_path, args.reco)
    
print 'Elapsed time:', datetime.datetime.now() - starttime

