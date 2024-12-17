import json
# from motify_functions import *

# un = {}
# un_en = {}
syn = {}
with open(f'all_syn.json', 'r') as f:
    datasyn = json.load(f)
    syn.update(datasyn)  

# fbl_items = {k: v for k, v in syn.items() if 'FBL' in k}
sorted_fbl_items = sorted(syn.items(), key=lambda item: item[1], reverse=True)
top_10_fbl = sorted_fbl_items[:10]

sorted_all_items = sorted(syn.items(), key=lambda item: item[1])
bottom_10_all = sorted_all_items[:10]

print("Top 10 FBL:", top_10_fbl)
print("Bottom 10 All:", bottom_10_all)
