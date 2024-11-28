from motify_functions import *
import sys
sys.path.append("../..")

file_dict = {}
directory = Path('cellcollective')
counter = 0
for file_path in directory.rglob('*'):
    if file_path.is_file():
        file_dict[counter] = str(file_path)
        counter += 1

functionsname = [fanout, fanin, cascade, mutualout, mutualin, bimutual, 
                 FFL, FBL, regulatingmutual, regulatedmutual, mutualcascade, semiclique, clique]

Apoptosis_un = {}
Apoptosis_un_en = {}
Apoptosis_syn = {}
for function in functionsname:
    un, un_en, syn = cal_motify(function, file_dict[0])
    Apoptosis_un.update(un)
    Apoptosis_un_en.update(un_en)
    Apoptosis_syn.update(syn)
x_values = list(Apoptosis_un.keys())
y_un = list(Apoptosis_un.values())
y_un_en = list(Apoptosis_un_en.values())
y_syn = list(Apoptosis_syn.values())
make_figure(x_values, y_un, y_un_en, y_syn, "Apoptosis_")

