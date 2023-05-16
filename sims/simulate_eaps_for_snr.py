"""
    Code should be run using:
        (1) Conda environment = n_prof_env


"""

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
import pickle

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
path2savemodels = os.path.join(workspace_dir,'Haynes2021_EAPs','model_data','bbp','snr_set')
path2figs = os.path.join(workspace_dir,'Haynes2021_EAPs','figures')


run_NEURON = True
pre_compiled = True
test_rest = False
save_df = True

vexc_factor = 1 # conversion factor for units

# for relative grid
pt3dRelativeToCellLocation = False # "Make grid as function from soma location (at origin) (default: True)"
invertedYCoord = False # " Make y-axis coordinate negative so they represent depth when visualized (0 at the top) (default: True)"




def valid_spike(Vm):

    if np.max(Vm)>0:
        return True
    else:
        return False



def simulate_eaps(model_num=0):

    # enforce integer
    if not type(model_num) == int:
        try:
            model_num = int(model_num)
        except ValueError:
            model_num = int(float(model_num))

    
    stim_comp = 'soma'
    secs = ['soma_0']
    ss_delay = 1000
    dur = 50
    pulse_dur = 10

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
    'NMLCL000651'

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
    'NMLCL000501'
    'NMLCL000499',
    'NMLCL000498',
    'NMLCL000497',
    'NMLCL000500'
    ]


    cell_labels = ['L5_TTPC1_%s'%(i+1) for i in range(5)]+  \
                  ['L5_TTPC2_%s'%(i+1) for i in range(5)]+   \
                  ['L5_UTPC_%s'%(i+1) for i in range(5)]+  \
                  ['L5_STPC_%s'%(i+1) for i in range(5)]+  \
                  ['L6_TPC_L4_%s'%(i+1) for i in range(5)]+\
                  ['L6_TPC_L1_%s'%(i+1) for i in range(5)]+\
                  ['L6_IPC_%s'%(i+1) for i in range(5)]+ \
                  ['L6_UTPC_%s'%(i+1) for i in range(5)]+ \
                  ['L6_BPC_%s'%(i+1) for i in range(5)]+ \
                  ['L4_SS_%s'%(i+1) for i in range(5)]+ \
                  ['L4_SP_%s'%(i+1) for i in range(5)]+ \
                  ['L23_NBC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_MC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_BTC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_DBC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_BPC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_LBC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_SBC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_ChC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_NGC_%s'%(i+1) for i in range(5)]
                

    

    pops1 = ['L5_PC' for _ in range(10)]
    pops2 = ['L6_PC' for _ in range(25)]
    pops3 = ['L4_PC' for _ in range(10)]
    pops4 = ['L23_IN' for _ in range(45)]
    pop_labels = pops1+pops2+pops3+pops4

    NMLDB_ID, cellLabel, popLabel = model_ids[model_num], cell_labels[model_num], pop_labels[model_num]

    from neuron import h
    # h.finitialize()

    example_neuron = utils.get_neuron_model_details(NMLDB_ID)

    file = example_neuron['model']['File_Name']
    model_dir = str(NMLDB_ID)

    if run_NEURON:
        fname = file[:-9] # Needed for mod in cellParams
        file =  fname + '.hoc'
        model_dir = model_dir + '-NEURON'

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
    rheobase = example_neuron['model']['Rheobase_High'] # choose _Low or _High
#     pulse_dur = 10
    pulse_threshold = example_neuron['model']['Threshold_Current_High']
    bias_curr = example_neuron['model']['Bias_Current'] # nA
    vrest = example_neuron['model']['Resting_Voltage']

    vinit = vrest
#     curr_amp = utils.get_short_square_current_amp(example_neuron)
    curr_amp = 3*rheobase





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
    y_shift = -500 # (microns)

    # get dimensions divisible by channel_spacing
    # (DEPRECATED)
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


    if run_NEURON:
        # set voltage initial conditions for all compartments
        for sec in importedCellParams['secs']:
            importedCellParams[sec]['vinit']=vinit # force to be same as in NEURON file


            
    target_pop=popLabel

    # drive spike to steady-state
    if popLabel in ['L23_PC']: bias_curr*=5 # bias current adjustment for L23_PC
    netParams.stimSourceParams['SSInput'] = {'type': 'IClamp',
                                                'del': 0, 'dur': ss_delay+dur,
                                                'amp': bias_curr}
    netParams.stimTargetParams['SSInput->'+cellLabel] = {'source': 'SSInput', 'sec':'soma_0', 'loc': 0.5,
                                                        'conds': {'pop':popLabel, 'cellList': [0]}}
    
    
    # apply injection at 3Xrheobase
    delay_onset = 0 # or 2
    netParams.stimSourceParams['Input'] = {'type': 'IClamp', 'del': ss_delay+delay_onset, 'dur': dur, 'amp': curr_amp}
    netParams.stimTargetParams['Input->'+cellLabel] = {'source': 'Input', 'sec':'soma_0', 'loc': 0.5,
                                                        'conds': {'pop':popLabel, 'cellList': [0]}}


    #### Define recording geometry
    r_pos = np.arange(10,125,5) # 21 r-values
    thetas = np.random.uniform(-np.pi,np.pi,size=10)
    
    all_xs = []
    all_zs = []
    
    for ri in r_pos:

        xi = np.multiply(ri,np.cos(thetas))
        zi = np.multiply(ri,np.sin(thetas))
        
        all_xs+=list(xi)
        all_zs+=list(zi)
    
    all_xs = np.array(all_xs)
    all_zs = np.array(all_zs)      
        
    all_ys = np.arange(-10,15,5) # can compute SNR/amplitude decay along the Y-axis

    rec_probes = []
    probe_list = []
    
    for xi, zi in zip(all_xs,all_zs):
        probe = [[xi,yi,zi] for yi in all_ys]
        probe_list.append(probe)

        rec_probes += probe

                
    rec_electrode = rec_probes
    num_channels, _ = np.shape(rec_electrode)
    
        
    # save probe arrangement
    file = '%s_recording_snr_probes.pkl'%cellLabel
    f = os.path.join(path2savemodels,file)
    
    with open(f, 'wb') as file:
        pickle.dump(probe_list, file,protocol=3)
        
        
            
            
    #######################
    #### Sim details ####
    #######################
#     desired_sampling_rate = 30e3 # (Hz) based on AIBS neuropixel data -- will just need to upsample
#     dt_sampling = 1./desired_sampling_rate

#     print('Dt optimal: %s | Dt sampling: %s'%(dt_optimal,dt_sampling))
#     raise Exception('Check')


    simConfig = specs.SimConfig()
    
    simConfig.pt3dRelativeToCellLocation = pt3dRelativeToCellLocation # "Make grid as function from soma location (at origin) (default: True)"
    simConfig.invertedYCoord = invertedYCoord # " Make y-axis coordinate negative so they represent depth when visualized (0 at the top) (default: True)"


    # global parameters
    simConfig.hParams= {'celsius': 34.0, # default is 6.3
                        'v_init':  vinit}
                        # 'v_init' : vrest}


    simConfig.duration = dur+ss_delay+1 # (ms)
    round_dt = float('%.3f'%dt_optimal) # near optimal
    simConfig.dt = round_dt
    simConfig.verbose = False
    simConfig.recordStep = 3*round_dt # multiple of optimal
    
    print('Sim dt: %s | Record dt: %s'%(simConfig.dt,simConfig.recordStep))


    simConfig.recordCells = ['all']  # which cells to record from
    simConfig.recordTraces = {'Vsoma':{'sec':'soma_0','loc':0.5,'var':'v'}}
    simConfig.recordStim = False
    simConfig.recordLFP = rec_electrode
    
    


    if save_df:
        (pops, cells, conns, stims, simData) = sim.createSimulateAnalyze(netParams = netParams, simConfig = simConfig,
                                output = True)

        print('Formatting data...')
        columns = ['Model_ID','t','vm','ve','x_bar','y_bar','z_bar',
                   'num_spikes','did_spike','first_spkt']
        waveforms_df = pd.DataFrame(columns=columns)

        # cast hoc objects to numpy.arrays
        try:
            did_spike = True
            num_spikes = len(np.array(simData['spkt']))
            spkt = np.array(simData['spkt'])[0]
            print('Yes spike!!')
        except IndexError:
            did_spike = False
            num_spikes = 0
            spkt = np.nan
            print('No spike!')

        t = np.array(simData['t'])
        Vm = np.array(simData['Vsoma']['cell_0'])
        Ve = np.array(simData['LFP'])

        temp_Ve = Ve.T


        # # visualization interval (defined per cell per condition)
        max_t = 1040 # never use 1050 because of simulation artefact at end

        temp_Ve = Ve.T
        rec_electrode = np.array(rec_electrode)

        for k in range(num_channels):
            
            # center point for electrode grid
            x_bar = rec_electrode[k,0]
            y_bar = rec_electrode[k,1]
            z_bar = rec_electrode[k,2]


            # grab current jitter Ve, reduce size and convert to uV
            ve_k = temp_Ve[k]
            
            this_interval = (t>=ss_delay)&(t<=max_t)
            
            try:
                reduced_ve_k = ve_k[this_interval] 
            except IndexError:
                this_interval = this_interval[:-1] # adjust by 1 timestep
                reduced_ve_k = ve_k[this_interval]
            
            try:
                reduced_t = t[this_interval]
                reduced_Vm = Vm[this_interval]
            except IndexError:
                this_interval = (t>=ss_delay)&(t<=max_t) # adjust it back
                reduced_t = t[this_interval]
                reduced_Vm = Vm[this_interval]
                

            # temp DataFrame
            df = pd.DataFrame(columns=columns)

            df['Model_ID'] = [NMLDB_ID]
            df['t'] = [reduced_t]
            df['vm'] = [reduced_Vm]
            df['ve'] = [reduced_ve_k]
            df['x_bar'] = [x_bar]
            df['y_bar'] = [y_bar]
            df['z_bar'] = [z_bar]
            df['first_spkt'] = [spkt]
            df['did_spike'] = [did_spike]
            df['num_spikes'] = [num_spikes]

            # join DataFrames
            join_frames = [waveforms_df,df]
            waveforms_df = pd.concat(join_frames,ignore_index=True)



        print('Saving data...')
        f = cellLabel + '_simulated_snr_eaps.pkl'
        filename = os.path.join(path2savemodels,f)  # Set file output name

        waveforms_df.to_pickle(filename,protocol=3)


        t = np.array(simData['t'])
        Vm = np.array(simData['Vsoma']['cell_0'])


        # test model
        result = valid_spike(Vm)

        if result:
            print('>>>>>>>>>>>>>>>>>>>>>>>>>>',NMLDB_ID,'>>',cellLabel,'>>',popLabel,'PASSED spike peak test!')
        else:
            print('>>>>>>>>>>>>>>>>>>>>>>>>>>',NMLDB_ID,'>>',cellLabel,'>>',popLabel,'FAILED spike peak test...')


    else:
         sim.createSimulateAnalyze(netParams = netParams, simConfig = simConfig,
                                output = False)

    



if __name__=='__main__':
    import sys


    simulate_eaps(*sys.argv[1:])
