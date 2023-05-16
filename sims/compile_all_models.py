# General modules
import os
import sys
import numpy as np
import matplotlib as plt
import pandas as pd
from pprint import pprint

# Simulation and analysis modules
from netpyne import sim, specs
import neuron

# relative import of NeuronProfiler
cwd = os.getcwd()
parent_dir = os.path.dirname(cwd)
workspace_dir = os.path.dirname(parent_dir)
sys.path.append(parent_dir)
sys.path.append(workspace_dir)
# print('Added\n --> %s \nto path.' %parent_dir)

import utils

model_set = 'test_bbp_models' # 'bbp' or 'aibs'

path2models = os.path.join(workspace_dir,'Haynes2021_EAPs',model_set)
path2savemodels = os.path.join(workspace_dir,'Haynes20201_EAPs','model_data','bbp')
path2figs = os.path.join(parent_dir,'figures')


run_NEURON = True
pre_compiled = False


model_ids = [
    # L5_TTPC1
    'NMLCL000687',
    'NMLCL000688',
    'NMLCL000689',
    'NMLCL000690',
    'NMLCL000691',

    # L5_TTPC2
    'NMLCL000695',
    'NMLCL000692', 
    'NMLCL000693', 
    'NMLCL000694', 
    'NMLCL000696', 

    # L5_UTPC
    'NMLCL000698',
    'NMLCL000699', 
    'NMLCL000701',
    'NMLCL000700',
    'NMLCL000697',

    # L5_STPC
    'NMLCL000684',
    'NMLCL000685',
    'NMLCL000682',
    'NMLCL000686',
    'NMLCL000683', 

    # L6_TPC_L4
    'NMLCL000672',
    'NMLCL000676',
    'NMLCL000673',
    'NMLCL000674', 
    'NMLCL000675',

    # L6_TPC_L1
    'NMLCL000669',
    'NMLCL000670',
    'NMLCL000667',
    'NMLCL000671',
    'NMLCL000668', 

    # L6_IPC
    'NMLCL000662',
    'NMLCL000666',
    'NMLCL000663',
    'NMLCL000664',
    'NMLCL000665',

    # L6_UTPC
    'NMLCL000680',
    'NMLCL000677',
    'NMLCL000681',
    'NMLCL000678',
    'NMLCL000679',

    # L6_BPC
    'NMLCL000659',
    'NMLCL000660', 
    'NMLCL000657',
    'NMLCL000661',
    'NMLCL000658',

    # L4_SS
    'NMLCL000652',
    'NMLCL000653',
    'NMLCL000654',
    'NMLCL000655',
    'NMLCL000656',

    # L4_SP
    'NMLCL000647',
    'NMLCL000648',
    'NMLCL000649',
    'NMLCL000650',
    'NMLCL000651',

    # L23_PC
    'NMLCL000637',
    'NMLCL000638',
    'NMLCL000639',
    'NMLCL000640',
    'NMLCL000641',

    # L23_NBC
    'NMLCL000494',
    'NMLCL000492',
    'NMLCL000496',
    'NMLCL000495', 
    'NMLCL000493', 

    # L23_MC
    'NMLCL000487',
    'NMLCL000491', 
    'NMLCL000488', 
    'NMLCL000489', 
    'NMLCL000490', 

    # L23_BTC
    'NMLCL000469',
    'NMLCL000487', 
    'NMLCL000491',
    'NMLCL000488',
    'NMLCL000490',

    # L23_DBC
    'NMLCL000480',
    'NMLCL000479',
    'NMLCL000477',
    'NMLCL000481', 
    'NMLCL000478',

    # L23_BP
    'NMLCL000462',
    'NMLCL000466',
    'NMLCL000463', 
    'NMLCL000464', 
    'NMLCL000465', 

    # L23_LBC
    'NMLCL000484',
    'NMLCL000485', 
    'NMLCL000482',
    'NMLCL000486',
    'NMLCL000483',

    # L23_SBC
    'NMLCL000505',
    'NMLCL000502',
    'NMLCL000506',
    'NMLCL000503',
    'NMLCL000504',

    # L23_ChC
    'NMLCL000473',
    'NMLCL000474', 
    'NMLCL000476',
    'NMLCL000472',
    'NMLCL000475', 

    # L23_NGC
    'NMLCL000501',
    'NMLCL000499',
    'NMLCL000498',
    'NMLCL000497',
    'NMLCL000500'
]


for NMLDB_ID in model_ids:
    print(NMLDB_ID)

    example_neuron = utils.get_neuron_model_details(NMLDB_ID)

    file = example_neuron['model']['File_Name']
    model_dir = str(NMLDB_ID)

    if run_NEURON:
        fname = file[:-9] # Needed for mod in cellParams
        file =  fname + '.hoc'
        model_dir = model_dir + '-NEURON'

    else:
        fname = file[:-9]

    model_path = os.path.join(path2models,model_dir)
    example_neuron_file = os.path.join(model_path,file)

    if not pre_compiled:
        print('Compiling model (.mod) files...')
        # Uncomment 'os.system()' on first run
        cmd_line_txt = '''
                      cd %s
                      nrnivmodl
                      ''' %model_path
        os.system(cmd_line_txt)
