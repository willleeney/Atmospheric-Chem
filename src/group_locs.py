
start_compounds = pickle.load(open('start_compounds.p', 'rb'))
compound_vps = start_compounds.copy()
base_id = list(np.zeros(len(compound_vps)))
locs = []

for idx, row in compound_vps.iterrows():
    cmp = row['Compound_Name']
    tloc =[1]
    for i in range(len(cmp)):
        if cmp[i] == 'C':# if you reach a new C 
            if (i+1) == len(cmp): # check for end of sequence 
                tloc += [1]
            elif cmp[i+1] == ')': # if add the end of a new chain 
                tloc += [1]
            elif cmp[i+1:i+5] == '(C)(' or  cmp[i+1] == '=' or cmp[i+1:i+3] == '(=':
                tloc += [0]
            elif cmp[i+1] == 'C': # if between two Cs
                if (i+2) == len(cmp):
                    tloc += [3]
                elif cmp[i+2] == '=' or cmp[i-1] == '=': # if one has a double bond 
                    tloc += [2]
                else:
                    tloc += [3] # if one doesnt   
            else:
                tloc += [2]
    locs += [tloc]

C_locs = []
for k in range(len(compound_vps)): 
    cp = compound_vps.loc[k, 'Compound_Name']
    c_add = [0]
    c_add += [idx+1 for idx,item in enumerate(cp) if item == 'C']
    C_locs += [c_add]

to_add = {'Compound_Name': list(compound_vps.iloc[:]['Compound_Name']), 'ID': base_id, 'Place': list(locs), 'Loc': C_locs}
to_merge = pd.DataFrame(to_add, columns = ['Compound_Name', 'ID', 'Place', 'Loc'])
compound_vps = pd.merge(compound_vps, to_merge)
pickle.dump(compound_vps, open("compound_vps.p", "wb"))
