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
path2savemodels = os.path.join(workspace_dir,'Haynes2021_EAPs','model_data','bbp')
path2figs = os.path.join(workspace_dir,'Haynes2021_EAPs','figures')


run_NEURON = True
pre_compiled = True
test_rest = False
save_df = True
sampling_dist = 'adjusted-uniform' # 'fixed-normal' or 'adjusted-lognormal' or 'adjusted-uniform'

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
    grid_spacing = 10 # 20 gives weird spread for BBP models, but technically is more accurate
    
    # total channels
#     channel_lb = 15
#     channel_ub = 16
#     channel_lb = 25
#     channel_ub = 25
    
    channel_lb = 32
    channel_ub = 32

    
    
    total_channels = channel_lb + channel_ub

#     model_ids = [
#         # Deep pyramids
#         'NMLCL000687', # works
#         'NMLCL000695',
#         'NMLCL000698',
#         'NMLCL000684',
#         'NMLCL000672',
#         'NMLCL000669',
#         'NMLCL000662',
#         'NMLCL000680',
#         'NMLCL000659',
#         # Superficial pyramids
#         'NMLCL000637',
#         'NMLCL000638',
#         'NMLCL000639',
#         'NMLCL000640',
#         'NMLCL000641',
#         # Interneurons
#         'NMLCL000494',
#         'NMLCL000487',
#         'NMLCL000469',
#         'NMLCL000480',
#         'NMLCL000462',
#         'NMLCL000484',
#         'NMLCL000505',
#         'NMLCL000473',
#         'NMLCL000501'
#     ]

#     cell_labels = [
#         # Deep pyramids
#         'L5_TTPC1',
#         'L5_TTPC2',
#         'L5_UTPC',
#         'L5_STPC',
#         'L6_TPC_L4',
#         'L6_TPC_L1',
#         'L6_IPC',
#         'L6_UTPC',
#         'L6_BPC',
#         # Superficial pyramids
#         'L23_PC1',
#         'L23_PC2',
#         'L23_PC3',
#         'L23_PC4',
#         'L23_PC5',
#         # Interneurons
#         'L23_NBC',
#         'L23_MC',
#         'L23_BTC',
#         'L23_DBC',
#         'L23_BP',
#         'L23_LBC',
#         'L23_SBC',
#         'L23_ChC',
#         'L23_NGC'
#     ]

#     pop_1 = ['L5_PC' for i in range(4)]
#     pop_2 = ['L6_PC' for i in range(5)]
#     pop_3 = ['L23_PC' for i in range(5)]
#     pop_4 = ['L23_IN' for i in range(9)]

#     pop_labels = pop_1+pop_2+pop_3+pop_4

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
                  ['L23_PC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_NBC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_MC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_BTC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_DBC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_BPC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_LBC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_SBC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_ChC_%s'%(i+1) for i in range(5)]+ \
                  ['L23_NGC_%s'%(i+1) for i in range(5)]

    
    
    pops1 = ['L5_PC' for _ in range(20)]
    pops2 = ['L6_PC' for _ in range(25)]
    pops3 = ['L4_PC' for _ in range(10)]
    pops4 = ['L23_PC' for _ in range(5)]
    pops5 = ['L23_IN' for _ in range(45)]
    pop_labels = pops1+pops2+pops3+pops4+pops5

    NMLDB_ID, cellLabel, popLabel = model_ids[model_num], cell_labels[model_num], pop_labels[model_num]
    
    print('Simulating %s'%cellLabel)

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

    

    #### Define recording electrodes geometry   
    num_probes = 100

    # sample radii from normal distribution
    print('Adjusting sample space: %s'%sampling_dist)
    if sampling_dist == 'fixed-normal':
        r_pos = np.random.normal(loc=20,scale=5,size=num_probes)
        
        min_radius = 10.
        r_pos = [r if r>min_radius else min_radius for r in r_pos] # force lower bound of 10.
        
    elif sampling_dist == 'adjusted-uniform':
        tempf = os.path.join(path2savemodels,'all_simulated_snr_dists.pkl')
        sampling_df = pd.read_pickle(tempf)
        

        df = sampling_df[sampling_df['cell_label']==cellLabel]
        empir_dist = df.empir_dist.iloc[0] # empir is a conservative estimate for most cases -- 20uV
        
        if empir_dist<=10:
            raise Exception("Model %s cannot be simulated -- model is dark neuron"%cellLabel) # ignore this model for adjusted dataset
        
        r_pos = np.random.uniform(low=10.,high=empir_dist,size=num_probes) # already includes upper and lower bound
        
        
    elif sampling_dist == 'adjusted-lognormal':
        tempf = os.path.join(path2savemodels,'all_simulated_snr_dists.pkl')
        sampling_df = pd.read_pickle(tempf)
        
        df = sampling_df[sampling_df['cell_label']==cellLabel]
        empir_dist = df.empir_dist.iloc[0] # empir is a conservative estimate for most cases -- 20uV
        
        if empir_dist<=10:
            raise Exception("Model %s cannot be simulated -- model dark neuron"%cellLabel) # ignore this model for adjusted dataset
        
        mu = np.log(empir_dist)
        sigma=0.3 # fixed spread
        r_pos = np.random.lognormal(mean=mu,sigma=sigma,size=num_probes)
        
        max_radius = theor_dist # force upper bound at noise limit
        min_r = np.min(r_pos)
        
        # shift really distant values within the range
        r_pos = [r if r<max_radius else np.random.uniform(min_r,max_radius-5,1)[0] for r in r_pos] # this is generally much less than 1% of units
                                                                                                    # indexing the random part to make sure proper data type
        
    else:
        raise Exception('Sampling distribution incorrectly specified.')
    

    # sample angle from uniform distribution
    theta_pos = np.random.uniform(-np.pi,np.pi,size=num_probes)

    all_xs = np.multiply(r_pos,np.cos(theta_pos))
    all_zs = np.multiply(r_pos,np.sin(theta_pos))

    # sample max channel height from normal distribution
#     all_ys = np.random.normal(0,1,size=num_probes) # sigma is small so misalignment may be 2-4 channels at most
#                                                    # upper and lower bound have 10+ channel buffer
#                                                    # needed because soma/AIS orientation differences alters max amplitude
#                                                    # channel, this is different across models and not known a priori

    all_ys = np.random.uniform(low=-0.5*grid_spacing,
                               high=0.5*grid_spacing,
                              size=num_probes)

    rec_probes = []
    probe_list = []
    
    lb = int(channel_lb*grid_spacing)
    ub = int(channel_ub*grid_spacing)
    
    # iterate over probe details
    for y_center, x_center, z_center in zip(all_ys,all_xs,all_zs):

        y_lower = np.arange(y_center-lb,y_center,grid_spacing)
        y_upper = np.arange(y_center,y_center+ub,grid_spacing)

        y_channels = list(y_lower) + list(y_upper)

        probe = [[x_center,yi,z_center] for yi in y_channels]
        probe_list.append(probe)

        rec_probes += probe
    
    rec_electrode = rec_probes
    
    num_channels, _ = np.shape(rec_electrode)
    
    # save probe arrangement
    file = '%s_opt_recording_probes.pkl'%cellLabel
    f = os.path.join(path2savemodels,file)
    
    with open(f, 'wb') as file:
        pickle.dump(probe_list, file,protocol=3)
        
        
            
            
    #######################
    ####  Sim details  ####
    #######################
        

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
        f = cellLabel + '_opt_simulated_eaps.pkl'
        filename = os.path.join(path2savemodels,f)  # Set file output name

        waveforms_df.to_pickle(filename,protocol=3)


#         print('TEST',simData)

        t = np.array(simData['t'])
        Vm = np.array(simData['Vsoma']['cell_0'])

#         f = cellLabel +'_vsoma_spike.png'
#         fig_f = os.path.join(path2figs,f)

#         fig = plt.figure(figsize=(8,8))
#         plt.plot(t,Vm)
#         print('Saving Vsoma plot')
#         fig.savefig(fig_f)

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
