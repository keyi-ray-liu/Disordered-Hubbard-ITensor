import sys
import json
from copy import deepcopy
import numpy as np
from itertools import product
import os
from shutil import copy2

catagories={
    #"sdcontact" : [True],
    #"sdcouple" : [2],
    "U" : [2.0, 3.0],
    #"TEcutoff": [ 1e-8 ],
    "L" : [128],
    "TEdim" : [128],
    "tswitch": [96],
    #"t_switch" : [10.0],
    #"fin" : [20.0]
    "mixed" : [True, False]
}

# sdpara = {
#     "to_chain_hop" : 0.0,
#     "internal_hop" : 1.0,
#     "source_offset" : 0.0,
#     "drain_offset" : 0.0,
#     "bulk_bias" : 0.0,
#     "init_bulk_bias" : [-100.0, 100.0],
#     "mixbasis" : False,
#     "NR" : 1,
#     "inelastic" : False,
#     "inelastic_para":{
#         "kcnt" : 0,
#         "deltak" : 0
#     },
#     "interacting" : False,
#     "sdcontact": True,
#     "Usd" : 0.0,
#     "sdcouple" : 2,
#     "state_override" : [2, 1]
# }

# transportpara = {
#     "t": 1,
#     "L": 2,
#     "N": [1, 1, 0, 0],
#     "sweepdim": 256,
#     "TEdim": 128,
#     "sweepcnt" : 50,
#     "TEcutoff" : 1e-9,
#     "krylovdim" : 10,
#     "CN" : 0,
#     "U" : 0.0,
#     "int_ee" : 0.0,
#     "int_ne" : 0.0,
#     "type" : "Fermion"
# }

dptpara = {
    "U" : 2.0,
    "L" : 4,
    "R" : 4,
    "tswitch" : 20.0,
    "timestep" : 0.25,
    "TEdim" : 32,
    "mixed" : False
}

keys = list(catagories.keys())
id = [ range(len(catagories[key])) for key in keys]



for prod in product(*id):

    # new_sd = deepcopy(sdpara)
    # new_transport = deepcopy(transportpara)
    new = deepcopy(dptpara)
    string =  '_'.join([ keys[i] + str(catagories[keys[i]][prod[i]]) for i in range(len(prod))]) 

    for key_i, val_i in enumerate(prod):

        changed = 0
        key = keys[key_i]
        option = catagories[key][val_i]

        # if key in new_sd:
        #     print(key, option)
        #     new_sd[key] = option
        #     changed = 1

        # if key in new_transport:
        #     print(key, option)
        #     new_transport[key] = option
        #     changed = 1

        if key in new:
            print(key, option)
            new[key] = option
            changed = 1

            if key == "L":
                new["R"] = option


        # if noneq:
        #     if key == 'Usd':

        #         new_sd['source_offset'] = option / 2 if option != 0 else 0.5
        #         new_sd['drain_offset'] = - option / 2 if option != 0 else -0.5

        #         changed = 1

        if changed == 0:
            raise(ValueError("category key '{}' does not match paras!".format(key)))

    if new["tswitch"] > new["L"]:
        continue 

    print(string)

    target_path = '/wrk/knl20/ITensor/newsrc/'
    #target_path = '/wrk/knl20/ITensor/ref/'
    #target_path = '/Users/knl20/Desktop/Code/TN/src/'

    if not os.path.exists(string):
        os.mkdir(string)

    for file in os.listdir(target_path):

        #print(target_path + file)
        if os.path.isfile(target_path + file):
            copy2( target_path+ file, string)

    with open(string +'/dptpara.json', 'w') as f:
        json.dump(new, f, indent=4)
    # with open(string +'/sdpara.json', 'w') as f:
    #     json.dump(new_sd, f, indent=4)

    # with open(string +'/transportpara.json', 'w') as f:
    #     json.dump(new_transport, f, indent=4)




