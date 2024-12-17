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

functionsname = [fanout, fanin, cascade, mutualout, mutualin, bimutual, 
                 FFL, FBL, regulatingmutual, regulatedmutual, mutualcascade, semiclique, clique]

all_un = {}
all_un_en = {}
all_syn = {}
keys = list(file_dict.keys())[76:77]  
#for key in file_dict.keys():
for key in keys:
    net_un = {}
    net_un_en = {}
    net_syn = {}
    for function in functionsname:
        un, un_en, syn = cal_motify(function, file_dict[key])
        net_un.update(un)
        net_un_en.update(un_en)
        net_syn.update(syn)
    make_un_bar_figure(net_un, f"{file_dict[key][15:-4]}")
    make_un_en_bar_figure(net_un_en, f"{file_dict[key][15:-4]}")
    make_syn_bar_figure(net_syn, f"{file_dict[key][15:-4]}")
    all_un.update(net_un)
    all_un_en.update(net_un_en)
    all_syn.update(net_syn)
    #print(f"{file_dict[key]}已完成")
with open('un77.json', 'w') as f:
    json.dump(all_un, f)
with open('un_en77.json', 'w') as f:
    json.dump(all_un_en, f)
with open('syn77.json', 'w') as f:
    json.dump(all_syn, f)
#make_un_bar_figure(all_un, "ALL_")
#make_un_en_bar_figure(all_un_en, "ALL_")
#make_syn_bar_figure(all_syn, "ALL_")





# all_un = {}
# all_un_en = {}
# all_syn = {}
# for key in file_dict.keys():
#     net_un = {}
#     net_un_en = {}
#     net_syn = {}
#     for function in functionsname:
#         un, un_en, syn = cal_motify(function, file_dict[key])
#         net_un.update(un)
#         net_un_en.update(un_en)
#         net_syn.update(syn)
#     all_un.update(net_un)
#     all_un_en.update(net_un_en)
#     all_syn.update(net_syn)
#     print(f"{file_dict[key]}已完成")
# make_un_bar_figure(all_un, "ALL_")
# make_un_en_bar_figure(all_un_en, "ALL_")
# make_syn_bar_figure(all_syn, "ALL_")


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
