import json
from motify_functions import *

un = {}
un_en = {}
syn = {}

lists = ["0_12", "14_17", "19_21", "24_27", "28_29", "30_33", "35", "41_43", "44", "45_51", "52", "53_58", "61_64", "66", "67_69", "70_75"]
for i in lists[0:1]:
    with open(f'un{i}.json', 'r') as f:
        dataun = json.load(f)
        un.update(dataun)
    with open(f'un_en{i}.json', 'r') as f:
        dataunen = json.load(f)
        un_en.update(dataunen)
    with open(f'syn{i}.json', 'r') as f:
        datasyn = json.load(f)
        syn.update(datasyn)  

make_un_bar_figure(un, "ALL_")
make_un_en_bar_figure(un_en, "ALL_")
make_syn_bar_figure(syn, "ALL_")