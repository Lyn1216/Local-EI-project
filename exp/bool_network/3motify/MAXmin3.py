import json
# from motify_functions import *

# un = {}
# un_en = {}
syn = {}

lists = ["0_12", "14_17", "19_21", "24_27", "28_29", "30_33", "35", "41_43", "44", "45_51", "52", "53_58", "61_64", "66", "67_69", "70_75"]
for i in lists[0:]:
    # with open(f'un{i}.json', 'r') as f:
    #     dataun = json.load(f)
    #     un.update(dataun)
    # with open(f'un_en{i}.json', 'r') as f:
    #     dataunen = json.load(f)
    #     un_en.update(dataunen)
    with open(f'syn{i}.json', 'r') as f:
        datasyn = json.load(f)
        syn.update(datasyn)  

fbl_items = {k: v for k, v in syn.items() if 'FBL' in k}
sorted_fbl_items = sorted(fbl_items.items(), key=lambda item: item[1], reverse=True)
top_10_fbl = sorted_fbl_items[:10]

sorted_all_items = sorted(syn.items(), key=lambda item: item[1])
bottom_10_all = sorted_all_items[:10]

print("Top 10 FBL:", top_10_fbl)
print("Bottom 10 All:", bottom_10_all)
