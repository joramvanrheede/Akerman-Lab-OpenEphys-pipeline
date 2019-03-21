#!/usr/bin/env python3

import ephysviz as viz
import numpy as np
from mltools import mlstudy as mls
from mltools import mlproc as mlp
import deepdish as dd
from matplotlib import pyplot as plt
import imp
import sys

chan_nr     = sys.argv[1]
data_dir    = sys.argv[2]

mlp.runProcess(
    'ephys.compute_templates', # processor name
    {"timeseries": data_dir + '/pre' + chan_nr + '.mda.prv',"firings":data_dir + '/firings' + chan_nr + '.mda.prv'}, # inputs
    {"templates_out": data_dir + '/templates' + chan_nr + '.npy'}, # outputs
    {"clip_size":100} # parameters
)

templates=mls.loadMdaFile(data_dir + '/templates' + chan_nr + '.npy')
viz.view_templates(templates)

mlp.runProcess(
    'ephys.compute_cross_correlograms',
    {"firings":data_dir + '/firings' + chan_nr + '.mda.prv'},
    {"correlograms_out":data_dir + '/autocorrelograms' + chan_nr + '.hdf5'},
    {"samplerate":30000,"max_dt_msec":50,"bin_size_msec":1,"mode":'autocorrelograms'}
)

X=dd.io.load(data_dir + '/autocorrelograms' + chan_nr + '.hdf5')
viz.view_cross_correlograms(X['correlograms'])

plt.show()
