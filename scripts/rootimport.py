'''Importing ROOT is no easy task (if you want argparse's help to work)'''

import sys

_argv = sys.argv
sys.argv = []
import ROOT
ROOT.TObject
sys.argv = _argv
del _argv

