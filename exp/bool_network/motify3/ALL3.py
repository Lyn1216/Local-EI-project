import json
import sys
sys.path.append("../..")
from func.motify_functions import *

un = {}
un_en = {}
syn = {}
dicts = {}

# lists = ["0_12", "14_17", "19_21", "24_27", "28_29", "30_33", "35", "41_43", "44", "45_51", "52", "53_58", "61_64", "66", "67_69", "70_75"]
# for i in lists[0:]:
#     with open(f'un{i}.json', 'r') as f:
#         dataun = json.load(f)
#         un.update(dataun)
#     with open(f'un_en{i}.json', 'r') as f:
#         dataunen = json.load(f)
#         un_en.update(dataunen)
#     with open(f'syn{i}.json', 'r') as f:
#         datasyn = json.load(f)
#         syn.update(datasyn)  
#     with open(f'dicts{i}.json', 'r') as f:
#         datadicts = json.load(f)
#         dicts.update(datadicts) 

# with open('all_un.json', 'w') as f:
#     json.dump(un, f)
# with open('all_un_en.json', 'w') as f:
#     json.dump(un_en, f)
# with open('all_syn.json', 'w') as f:
#     json.dump(syn, f)
# with open('all_dicts.json', 'w') as f:
#     json.dump(dicts, f)

with open('all_un.json', 'r') as f:
    un = json.load(f)
with open('all_un_en.json', 'r') as f:
    un_en = json.load(f)
with open('all_syn.json', 'r') as f:
    syn = json.load(f)
with open('all_dicts.json', 'r') as f:
    dicts = json.load(f)
    
motifyname = ["fanout", "fanin", "cascade", "mutualout", "mutualin", "bimutual", "FFL", "FBL", "regulatingmutual", "regulatedmutual", "mutualcascade", "semiclique", "clique"]
make_un_bar_figure(un, "ALL_", motifyname)
make_un_en_bar_figure(un_en, "ALL_", motifyname)
make_syn_bar_figure(syn, "ALL_", motifyname)
make_combined_bar_figure(un, un_en, syn, 'ALL_combined_figure', motifyname)