# General modules
import os
import matplotlib as mpl
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
import matplotlib.pyplot as plt

import sys
import numpy as np
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

runModelsFrom = 'aibs'

if runModelsFrom == 'bbp':
    model_set = 'test_bbp_models' # 'bbp' or 'aibs'
    
    model_ids = [
        # Deep pyramids
        'NMLCL000687', 
        'NMLCL000695',
        'NMLCL000698',
        'NMLCL000684',
        'NMLCL000672',
        'NMLCL000669',
        'NMLCL000662',
        'NMLCL000680',
        'NMLCL000659',
        # Superficial pyramids
        'NMLCL000637',
        'NMLCL000638',
        'NMLCL000639',
        'NMLCL000640',
        'NMLCL000641',
        # Interneurons
        'NMLCL000494',
        'NMLCL000487',
        'NMLCL000469',
        'NMLCL000480',
        'NMLCL000462',
        'NMLCL000484',
        'NMLCL000505',
        'NMLCL000473',
        'NMLCL000501'
    ]

    cell_labels = [
        # Deep pyramids
        'L5_TTPC1',
        'L5_TTPC2',
        'L5_UTPC',
        'L5_STPC',
        'L6_TPC_L4',
        'L6_TPC_L1',
        'L6_IPC',
        'L6_UTPC',
        'L6_BPC',
        # Superficial pyramids
        'L23_PC1',
        'L23_PC2',
        'L23_PC3',
        'L23_PC4',
        'L23_PC5',
        # Interneurons
        'L23_NBC',
        'L23_MC',
        'L23_BTC',
        'L23_DBC',
        'L23_BP',
        'L23_LBC',
        'L23_SBC',
        'L23_ChC',
        'L23_NGC'
    ]

    pop_1 = ['L5_PC' for i in range(4)]
    pop_2 = ['L6_PC' for i in range(5)]
    pop_3 = ['L23_PC' for i in range(5)]
    pop_4 = ['L23_IN' for i in range(9)]

    pop_labels = pop_1+pop_2+pop_3+pop_4
    
    path2models = os.path.join(workspace_dir,'Haynes2021_EAPs',model_set)
    path2savemodels = os.path.join(parent_dir,'model_data',model_set)
    path2figs = os.path.join(parent_dir,'figures')


if runModelsFrom == 'aibs':
    model_set = 'test_aibs_models'
    fname = 'aibs_model_types.csv'
    
    path2models = os.path.join(workspace_dir,'Haynes2021_EAPs',model_set)
    path2savemodels = os.path.join(parent_dir,'model_data',runModelsFrom)
    path2figs = os.path.join(parent_dir,'figures')
    
    f = os.path.join(path2models,fname)
    cell_tags_df = pd.read_csv(f)
    
    model_ids = cell_tags_df.id.tolist()
    cre_lines = cell_tags_df.cre_lines.tolist()
    model_types = cell_tags_df.model_type.tolist()
    
    cell_labels = []   
    pop_labels = []
    for cre, model_type in zip(cre_lines,model_types):
        tag1 = cre.split('-')[0]
        
        if model_type == 'perisomatic':
            tag2 = '_ps'
        else: #'all-active'
            tag2 = '_aa'
        
        cell_labels.append(tag1+tag2) # cre line + type of model
        pop_labels.append(tag1) # just cre line

    

    





runNEURON = False
runSONATA = True

pre_compiled = True
test_rest = False

def valid_spike(Vm):

    if np.max(Vm)>0:
        return True
    else:
        return False



def test_model(model_num=0):

    # enforce integer
    if not type(model_num) == int:
        try:
            model_num = int(model_num)
        except ValueError:
            model_num = int(float(model_num))


    stim_comp = 'soma'
    secs = ['soma_0']
    ss_delay = 0
    dur = 1000





    NMLDB_ID, cellLabel, popLabel = model_ids[model_num], cell_labels[model_num], pop_labels[model_num]

    from neuron import h
    # h.finitialize()

    example_neuron = utils.get_neuron_model_details(NMLDB_ID)

    file = example_neuron['model']['File_Name']
    model_dir = str(NMLDB_ID)

    if runNEURON:
        fname = file[:-9] # Needed for mod in cellParams
        file =  fname + '.hoc'
        model_dir = model_dir + '-NEURON'
        
    if runSONTATA:
        pass
    
    else:
        fname = file[:-9]
        file = fname+'.cell.nml'

    model_path = os.path.join(path2models,model_dir)
    example_neuron_file = os.path.join(model_path,file)

    if not pre_compiled:
        # print('Compiling model (.mod) files...')
        # Uncomment 'os.system()' on first run
        cmd_line_txt = '''
                      cd %s
                      nrnivmodl
                      ''' %model_path
        os.system(cmd_line_txt)

    # print('Loading mechanisms...')
    neuron.load_mechanisms(model_path)





    # Model run details
    dt_optimal = example_neuron['model']['Optimal_DT']
#     rheobase = example_neuron['model']['Rheobase_High'] # choose _Low or _High
#     pulse_dur = 10
    rheobase = example_neuron['model']['Threshold_current_High']
    pulse_dur = 3
    bias_curr = example_neuron['model']['Bias_Current'] # nA
    vrest = example_neuron['model']['Resting_Voltage']

    vinit = vrest
    curr_amp = utils.get_short_square_current_amp(example_neuron)






    netParams = specs.NetParams()


    ####################
    # Geometry
    # --------
    # (Cartesian axes)
    #          y    z
    #          ^  ^
    #          | /
    #   x <--- o --
    #         /|
    ####################

    netParams.propVelocity = 100.0 # propagation velocity (um/ms)
    netParams.probLengthConst = 150.0 # length constant for conn probability (um)

    # retrieve model height (to place pseudo-electrode points)
    example_neuron_morpho = utils.get_neuron_model_morphometrics(NMLDB_ID)
    for metric_dict in example_neuron_morpho:
        if metric_dict['Metric_ID'] == 'Height':
            model_height = metric_dict['Maximum']
        if metric_dict['Metric_ID'] == 'Width':
            model_width = metric_dict['Maximum']
        if metric_dict['Metric_ID'] == 'Depth':
            model_depth = metric_dict['Maximum']


    # define recording electrode geometry (get more precise)
    channel_spacing = 50. # (microns)
    buffer_dim = 50. # (microns)

    # TODO: Should be defined based on y minimum with soma placed at origin
    y_shift = -225 # (microns)

    # get dimensions divisible by channel_spacing
    x_dim =  channel_spacing*np.floor((model_width+buffer_dim)/channel_spacing)+channel_spacing
    y_dim =  channel_spacing*np.floor((model_height+buffer_dim)/channel_spacing)+channel_spacing
    z_dim =  channel_spacing*np.floor((model_depth+buffer_dim)/channel_spacing)+channel_spacing

    # make x-z plane square
    if x_dim>z_dim:
        z_dim = x_dim
    else:
        x_dim = z_dim

    netParams.sizeX = x_dim # x-dimension (horizontal length) size in um
    netParams.sizeY = y_dim # y-dimension (vertical height or cortical depth) size in um
    netParams.sizeZ = z_dim # z-dimension (horizontal length) size in um



    netParams.popParams[popLabel] = {'cellModel' : 'Mark2015',
                                     'cellType' : fname,
                                     'numCells' : 1}

    # import cellParams
    loadCellParams = NMLDB_ID

    importedCellParams = netParams.importCellParams(label=cellLabel,
                                somaAtOrigin=True,
                                conds={'cellType': fname, 'cellModel': 'Mark2015'},
                                fileName=example_neuron_file, cellName = fname)


    if runNEURON:
        # set voltage initial conditions for all compartments
        for sec in importedCellParams['secs']:
            importedCellParams[sec]['vinit']=vinit # force to be same as in NEURON file

    if not test_rest:
        # apply bias adjusted current injection at rheobase
        netParams.stimSourceParams['Input'] = {'type': 'IClamp', 'del': 10, 'dur': pulse_dur, 'amp': curr_amp}
        netParams.stimTargetParams['Input->'+cellLabel] = {'source': 'Input', 'sec':'soma_0', 'loc': 0.5,
                                                            'conds': {'pop':popLabel, 'cellList': [0]}}
    else:
        # apply bias adjusted current injection
        if popLabel in ['L23_PC']: bias_curr*=5 # bias current adjustment for L23_PC
        netParams.stimSourceParams['Input'] = {'type': 'IClamp', 'del': ss_delay, 'dur': dur, 'amp': bias_curr}
        netParams.stimTargetParams['Input->'+cellLabel] = {'source': 'Input', 'sec':'soma_0', 'loc': 0.5,
                                                            'conds': {'pop':popLabel, 'cellList': [0]}}

    #
    # else:
    # apply default current to achieve steady state
    # netParams.stimSourceParams['SS_Input'] = {'type': 'IClamp', 'del': 0, 'dur': dur+ss_delay, 'amp': 0}
    # netParams.stimSourceParams['SS_Input'] = {'type': 'IClamp', 'del': 0, 'dur': dur+ss_delay, 'amp': bias_curr}
    # netParams.stimTargetParams['SS_Input->'+cellLabel] = {'source': 'SS_Input', 'sec': 'soma_0', 'loc': 0.5,
    #                                                     'conds': {'pop':popLabel, 'cellList': [0]}}



    # This model is careful to ensure the clamp current is properly computed relative to the membrane voltage
    # netParams.stimSourceParams['SS_Voltage'] = {'type': 'SEClamp', 'dur1':0}
    # netParams.stimTargetParams['SS_Voltage->'+cellLabel] = {'source': 'SS_Voltage', 'sec': 'soma_0', 'loc': 0.5,
    #                                                     'conds': {'pop':popLabel, 'cellList': [0]}}


    # CHECK
    # f = open('test_secs_valid.txt','w')
    # f.write(str(netParams.cellParams[cellLabel]['secs']))
    # f.close()
    #
    # f = open('test_popparams_valid.txt','w')
    # f.write(str(netParams.popParams))
    # f.close()
    #
    # f = open('test_stimtargetparams_valid.txt','w')
    # f.write(str(netParams.stimTargetParams))
    # f.close()
    # raise Exception

    #######################
    #### Sim details ####
    #######################
    # ds_factor = 10 # downsampling factor

    simConfig = specs.SimConfig()

    # global parameters
    simConfig.hParams= {'celsius': 34.0, # default is 6.3
                        'v_init':  vinit}
                        # 'v_init' : vrest}


    simConfig.duration = dur+ss_delay # (ms)
    simConfig.dt = dt_optimal
    simConfig.verbose = True
    # simConfig.recordStep = ds_factor * dt_optimal
    simConfig.recordStep = 0.1


    simConfig.recordCells = ['all']  # which cells to record from
    simConfig.recordTraces = {'Vsoma':{'sec':'soma_0','loc':0.5,'var':'v'}}
    simConfig.recordStim = False


    (pops, cells, conns, stims, simData) = sim.createSimulateAnalyze(netParams = netParams, simConfig = simConfig,
                            output = True)

    print('TEST',simData)

    t = np.array(simData['t'])
    Vm = np.array(simData['Vsoma']['cell_0'])

    if test_rest: f = cellLabel +'_vsoma_rest_mem_test.png'
    else: f = cellLabel +'_vsoma_spike_test.png'
    fig_f = os.path.join(path2figs,f)

    fig = plt.figure(figsize=(8,8))
    plt.plot(t,Vm)
    print('Saving Vsoma plot')
    fig.savefig(fig_f)

    # test model
    result = valid_spike(Vm)

    if result:
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>',NMLDB_ID,'>>',cellLabel,'>>',popLabel,'PASSED spike peak test!')
    else:
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>',NMLDB_ID,'>>',cellLabel,'>>',popLabel,'FAILED spike peak test...')


    if test_rest:
        print('Steady State at %s vs final Vm at %s'%(vrest,Vm[-1]))






if __name__=='__main__':
    import sys


    test_model(*sys.argv[1:])
