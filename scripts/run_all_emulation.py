#!/usr/bin/env python
'''Run reconstruction emulation for all (selected) event samples.

Adapted from run_all_emulations.py by D. Kaleko.

Usage:

    python ./run_all_emulations.py [input_file_base] [output_file_base]

N.B. This calls out to another script with the shell.

Author: A. Mastbaum <mastbaum@uchicago.edu>, 2016/11/16
'''

import argparse
import os
import shlex
import subprocess

suffix = '_recoemmaster_'

parser = argparse.ArgumentParser(description='Run reconstruction emulation')
parser.add_argument('out_path', default='.', help='Output file path')
parser.add_argument('in_files', nargs='+', default='.', help='Input files')
args = parser.parse_args()

dev_path = os.getenv('LARLITE_USERDEVDIR')
module_path = os.path.join(dev_path, 'LArLiteApp', 'RecoEmulatorApp', 'mac')

script = os.path.join(module_path, 'run_reco_emulator.py')
config = os.path.join(module_path, 'recoemu_master.fcl')

for f in args.in_files:
    outfile = os.path.splitext(os.path.basename(f))[0] + suffix
    outfile = os.path.join(args.out_path, outfile)

    cmd = 'python %s %s %s %s' % (script, config, f, outfile)
    print 'Emulation command: %s' % cmd
    subprocess.call(shlex.split(cmd))

