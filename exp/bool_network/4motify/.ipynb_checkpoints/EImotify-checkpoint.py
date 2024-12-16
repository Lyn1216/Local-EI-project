from motify_functions import *
import sys
sys.path.append("../..")
import json

file_dict = {}
directory = Path('cellcollective')
counter = 0
for file_path in directory.rglob('*'):
    if file_path.is_file():
        file_dict[counter] = str(file_path)
        counter += 1

functionsname = [bifan, twoFFL, fourone, fourtwo, fiveone, fivetwo, fivethree, sixone, sixtwo, sixthree, sixfour, 
                 sevenone, seventwo, seventhree, eightone, eighttwo, nine, ten, eleven, twelve]

all_un = {}
all_un_en = {}
all_syn = {}
all_dicts = {}
keys = list(file_dict.keys())[24:28]  
for key in keys:
    net_un = {}
    net_un_en = {}
    net_syn = {}
    net_dicts = {}
    for function in functionsname:
        un, un_en, syn, dicts = cal_motify(function, file_dict[key])
        net_un.update(un)
        net_un_en.update(un_en)
        net_syn.update(syn)
        net_dicts.update(dicts)
    make_un_bar_figure(net_un, f"{file_dict[key][15:-4]}")
    make_un_en_bar_figure(net_un_en, f"{file_dict[key][15:-4]}")
    make_syn_bar_figure(net_syn, f"{file_dict[key][15:-4]}")
    all_un.update(net_un)
    all_un_en.update(net_un_en)
    all_syn.update(net_syn)
    all_dicts.update(net_dicts)
with open('un24_27.json', 'w') as f:
    json.dump(all_un, f)
with open('un_en24_27.json', 'w') as f:
    json.dump(all_un_en, f)
with open('syn24_27.json', 'w') as f:
    json.dump(all_syn, f)
all_dicts = {str(key): value for key, value in all_dicts.items()}
with open('dicts24_27.json', 'w') as f:
    json.dump(all_dicts, f)

# Arabidopsis_un = {}
# Arabidopsis_un_en = {}
# Arabidopsis_syn = {}
# for function in functionsname:
#     un, un_en, syn = cal_motify(function, file_dict[1])
#     Arabidopsis_un.update(un)
#     Arabidopsis_un_en.update(un_en)
#     Arabidopsis_syn.update(syn)

# make_un_bar_figure(Arabidopsis_un, "Arabidopsis_", file_dict[1])
# make_un_en_bar_figure(Arabidopsis_un_en, "Arabidopsis_", file_dict[1])
# make_syn_bar_figure(Arabidopsis_syn, "Arabidopsis_", file_dict[1])
