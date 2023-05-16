import os 
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

these_models = range(len(model_ids))
print('Simulating %s Models'%len(model_ids))

for model_num in these_models:
    
    cmd_line_txt = """
               python simulate_eaps.py %s
               """%model_num
    
    os.system(cmd_line_txt)
