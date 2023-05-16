# general modules
from ast import literal_eval
import pandas as pd
import numpy as np
import requests
import xml.etree.ElementTree as ET
import os
from os.path import join, isfile

def confirm_cortical_model(neurolex_term):
    '''
        Used to NeuroLex ID containing the brain area. 'cort' used for neocortex or other similar names.
        ----
        PARAMETERS : ()
        OUTPUT : (bool) Wheter a cortical model was found or not
    '''
    try:
        return 'cort' in neurolex_term
    except TypeError: # if neurolex_term is None, these should be manually validated
        return False

def get_model(nmldb_id,save=False,path_to_save=None):
    '''
        NeuroML-DB API query for model
        TODO: Check if file is store already somewhere and load that instead of calling API.
        ----
        PARAMETER:
            - nmldb_id : (str) Contains NML ID
            - save : (bool)
            - path_to_save : (str)
        OUTPUT:
            - mldb_model_response : (str) Contains all model specifications
    '''


    nmldb_xml_url = 'https://neuroml-db.org/render_xml_file?modelID='

    model_xml_url = nmldb_xml_url + nmldb_id

    xml_name = nmldb_id+'.xml'

    filename = join(path_to_save,xml_name)

    # If file exists, load and return
    if os.path.isfile(filename):
        return None # replace with uploading the model


    model_xml_response = requests.get(model_xml_url)


    # save file
    if save:

        if path_to_save:
            filename = os.path.join(path_to_save,xml_name)
        else:
            filename = xml_name

        with open(filename, 'wb') as file:
            file.write(model_xml_response.content)

    return model_xml_response.text





def get_neuron_model_details(nmldb_id):
    '''
        NeuroML-DB API query for model details
        TODO: Check if file is store already somewhere and load that instead of calling API.
        ----
        PARAMETER:
            - nmldb_id : (str) Contains NML ID
        OUTPUT:
            - mldb_model_response : (str) Contains all model details in JSON format
    '''

    # collects all model details from database
    if nmldb_id in ['all']:
        nmldb_url = 'http://neuroml-db.org/api/models'

        nmldb_model_response = requests.get(nmldb_url)

    # collects single model details
    else:
        nmldb_url = 'http://neuroml-db.org/api/model?id='

        neuron_url = nmldb_url + nmldb_id

        nmldb_model_response = requests.get(neuron_url)

    return nmldb_model_response.json()


def get_channel_model_details(nmldb_id):

    if 'CL' in nmldb_id:
        raise Exception('NMLDB ID corresponds to a cell model!')

    elif nmldb_id in ['all']:
        nmldb_url = 'http://neuroml-db.org/api/models'

        nmldb_model_response = requests.get(nmldb_url)

    else:
        nmldb_url = 'http://neuroml-db.org/api/model?id='

        neuron_url = nmldb_url + nmldb_id

        nmldb_model_response = requests.get(neuron_url)

    return nmldb_model_response.json()




def get_neuron_model_morphometrics(nmldb_id):

    nmldb_url = 'http://neuroml-db.org/api/morphometrics?id='

    neuron_url = nmldb_url + nmldb_id

    nmldb_morpho_response = requests.get(neuron_url)

    return nmldb_morpho_response.json()


def get_short_square_current_amp(example_neuron):

    protocol_curr_list = []
    for wave_dict in example_neuron['waveform_list']:


        if wave_dict['Protocol_ID'] == 'SHORT_SQUARE':
            protocol_curr_list.append(wave_dict)

    return float(protocol_curr_list[-1]['Waveform_Label'].split(' ')[0])



def get_model_waveform(waveform_id):

        waveform_url = 'http://neuroml-db.org/api/waveform?id=%s' %waveform_id

        nmldb_waveform_response = requests.get(waveform_url)

        return nmldb_waveform_response.json()


def import_model(filename,path_to_models=None):

    # first check if model file exists
    model_file = os.path.join(path_to_models,filename)
    print(model_file)

    if os.path.isfile(model_file):

        model = ET.parse(model_file)
        return model

    else:
        print('No model with NMLDB ID: %s exists!' %nmldb_id)

def query_model_by_keyword(keyword):

    nmldb_url = 'http://neuroml-db.org/api/search?q='

    keyword_url = nmldb_url + keyword

    nmldb_keyword_response = requests.get(keyword_url)


    return nmldb_keyword_response.json()
