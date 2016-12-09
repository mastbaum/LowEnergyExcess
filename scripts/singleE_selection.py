#!/usr/bin/env python

'''Perform single-electron event selection with various filters

Adapted from singleE_*_selection.py, by D. Kaleko

Usage:

    python singleE_selection.py [mc|reco] [sample] [input files] [output path]

A. Mastbaum <mastbaum@uchicago.edu>, 2016/11/16
'''

import argparse
import os
import sys

from ROOT import larlite as fmwk
from ROOT import ertool
from singleE_config import GetERSelectionInstance

filters = {
    'nue':             (fmwk.MC_CCnue_Filter,  'beamNuE'),
    'numu':            (fmwk.MC_CCnumu_Filter, 'beamNuMu'),
    'nc':              (fmwk.MC_NC_Filter,     'beamNC'),
    'cosmic':          (fmwk.MC_cosmic_Filter, 'cosmicShowers'),
    'dirt':            (fmwk.MC_dirt_Filter,   'dirt'),
    'LEE':             (fmwk.MC_LEE_Filter,    'LEETree'),
    'cosmicoutoftime': (None,                  'cosmicOutOfTime'),
}


def run(sample, infiles, outpath, reco=False):
    proc = fmwk.ana_processor()
    proc.enable_filter(True)

    for f in infiles:
        proc.add_input_file(f)

    # Specify IO mode
    proc.set_io_mode(fmwk.storage_manager.kBOTH)

    # Specify output ROOT file name
    filename = 'singleE_%s_selection_%s' % (sample, 'reco' if reco else 'mc')
    outfile = os.path.join(outpath, filename)

    proc.set_ana_output_file(outfile + '.root')
    proc.set_output_file(outfile + '_larlite_out.root')

    filter_type, treename = filters[sample]

    if filter_type is not None:
        event_filter = filter_type()
        proc.add_process(event_filter)

    ana = ertool.ERAnaLowEnergyExcess()
    ana.SetTreeName(treename)

    if sample == 'LEE':
        ana.SetLEESampleMode(True)

        # Set number of LEE events for the MC, after filtering
        ana.SetLEENEvents(13302)

        lee_file = \
            os.path.join(os.environ['LARLITE_USERDEVDIR'], 'LowEnergyExcess',
                         'LEEReweight', 'source', 'LEE_Reweight_plots.root')
        ana.SetLEEFilename(lee_file)

        # Currently using the mcc6 input histogram even though it's not quite
        # right, because I can't get the mcc7 input histogram to work
        # correctly with these low statistics - DK
        ana.SetLEECorrHistName('initial_evis_uz_corr')

    anaunit = GetERSelectionInstance()
    anaunit._mgr.ClearCfgFile()

    cfg_file = 'ertool_default_emulated.cfg' if reco else 'ertool_default.cfg'
    cfg_path = os.path.join(os.environ['LARLITE_USERDEVDIR'],
                            'SelectionTool', 'ERTool', 'dat', cfg_file)

    anaunit._mgr.AddCfgFile(cfg_path)

    if reco:
        anaunit.SetShowerProducer(False, 'recoemu')
        anaunit.SetTrackProducer(False, 'recoemu')

    anaunit._mgr.AddAna(ana)

    proc.add_process(anaunit)

    proc.run()


if __name__ == '__main__':
    desc = 'Single electron sample selections for the LEE analysis'
    parser = argparse.ArgumentParser(description=desc)
    parser.add_argument('-r', '--reco', action='store_true',
                        help='Use reconstructed quantities')
    parser.add_argument('sample', choices=filters.keys(),
                        help='Sample name to select')
    parser.add_argument('outpath', help='Output file path')
    parser.add_argument('inputs', nargs='+', help='Input files')
    args = parser.parse_args()

    run(args.sample, args.inputs, args.outpath, reco=args.reco)

