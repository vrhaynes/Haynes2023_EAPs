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
