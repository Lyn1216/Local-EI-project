import json
from motify_functions import *

un = {}
un_en = {}
syn = {}
dicts = {}

# lists = ["0", "1_3", "4_7", "8_12", "14_17", "19_21", "24_27", "28_29", "30_33", "35", "41_58", "61_64", "66_75"]
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

# make_un_bar_figure(un, "ALL_")
# make_un_en_bar_figure(un_en, "ALL_")
# make_syn_bar_figure(syn, "ALL_")
make_combined_bar_figure(un, un_en, syn, "ALL_combine_")

