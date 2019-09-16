################################## ADDING PARENT COMPOUND ID ######################################################
compound_vps = pickle.load(open('save_data/compound_vps.p', 'rb'))
start_compounds = pickle.load(open('save_data/start_compounds.p', 'rb'))
functional_groups = ['O', 'OO', 'C=O', 'C(=O)O', 'C#N', 'N', 'Cl', 'Br', 'I'] # each of the initial fun groups
f_group_start = ['O', 'OO', 'C(=O)', 'OC(=O)', 'N#C', 'N', 'Cl', 'Br', 'I']
ind = []
new_cmps = []
IDS = []
fun_addition = compound_vps.copy()

for i in tqdm(range(len(compound_vps))): # all framework compounds
    temp_cmp = compound_vps.at[i, 'Compound_Name']
    for j in range(len(temp_cmp)): # each position in compound
        if (j == 0):
            for k in range(len(functional_groups)): # each functional group
                cmp = '{}{}'.format(f_group_start[k], temp_cmp[:]) # N bond is different because is a triple bond 
                ind.append(int(i))
                new_cmps.append(cmp)
                IDS.append(k+1)
        elif(j == len(temp_cmp)-1):
            for k in range(len(functional_groups)): # each functional group
                cmp = '{}({}){}'.format(temp_cmp[:j+1], functional_groups[k], temp_cmp[j+1:])
                ind.append(int(i))
                new_cmps.append(cmp)
                IDS.append(k+1)
        elif temp_cmp[j+1] == ')':
            for k in range(len(functional_groups)): # each functional group
                cmp = '{}({}){}'.format(temp_cmp[:j+1], functional_groups[k], temp_cmp[j+1:])
                ind.append(int(i))
                new_cmps.append(cmp)
                IDS.append(k+1)
                
                
to_add = {'Compound_Name': new_cmps, 'Parent': ind, 'ID': IDS}
merge = pd.DataFrame(to_add, columns = ['Compound_Name', 'Parent', 'ID'])
fun_addition = pd.merge(fun_addition, merge, how='outer')#, left_on = ['Colour'], right_on = ['Vapour_Pressure'])

########### ASSIGNS INITIAL PARENT ID ###########
for i in range(1175):
    fun_addition.at[i, 'Parent'] = int(i)
#
############# CHANGES COLOUR FOR NEW FUN GROUPS ##########
#for index, row in fun_addition.iterrows():
#    temp = row['Compound_Name']
#    col = row['ID']
#    if col == 1: 
#        if 'COO' in temp: 
#            fun_addition.at[index, 'ID'] = 10
#    if col == 2: # this one doesnt work 
#        if 'OOCOO' in temp: 
#            fun_addition.at[index, 'ID'] = 11
#    if col == 4: 
#        if '(OCOO)' in temp:
#            fun_addition.at[index, 'ID'] = 12

 ######### ADDS PARENT MOLECULE FOR LONGER FUN GROUPS ####### this could be improved
example1 = fun_addition[fun_addition['ID'] == 1].copy()
example2 = fun_addition[fun_addition['ID'] == 2].copy()
example3 = fun_addition[fun_addition['ID'] == 4].copy()
#example = pd.merge(example, example1, how = 'outer')
for index, row in example1.iterrows():
    temp = row['Compound_Name']
    for j in range(len(temp)):
        if temp[j] == 'O':
            example1.at[index, 'Compound_Name'] = '{}{}{}'.format(temp[:j+1], 'COO', temp[j+1:])
            #example1.at[index, 'Colour'] = row['Parent']
            break 
            
for index, row in example2.iterrows():
    temp = row['Compound_Name']
    for j in range(len(temp)):
        if temp[j] == 'O':
            example2.at[index, 'Compound_Name'] = '{}{}{}'.format(temp[:j], 'OOC', temp[j:])
            #example2.at[index, 'Colour'] = row['Parent']
            break
                    
for index, row in example3.iterrows():
    temp = row['Compound_Name']
    for j in range(len(temp)):
        if temp[j:j+2] == '=O':
            example3.at[index, 'Compound_Name'] = '{}{}{}'.format(temp[:j+3], '(OCOO)', temp[j+4:])
            #example3.at[index, 'Colour'] = row['Parent']
            break 
        
fun_addition = fun_addition.append(example1, ignore_index=True)
fun_addition = fun_addition.append(example2, ignore_index=True)
fun_addition = fun_addition.append(example3, ignore_index=True)


############ CHANGES COLOUR FOR NEW FUN GROUPS ##########
for index, row in fun_addition.iterrows():
    temp = row['Compound_Name']
    col = row['ID']
    if col == 1: 
        if 'COO' in temp: 
            fun_addition.at[index, 'ID'] = 10
    if col == 2: # this one doesnt work 
        if 'OOCOO' in temp: 
            fun_addition.at[index, 'ID'] = 11
    if col == 4: 
        if '(OCOO)' in temp:
            fun_addition.at[index, 'ID'] = 12


#for adding an oxygen between two carbons in each and every valid place, then all places

new_cmps = []
parents = []
ids = []
for index, row in start_compounds.iterrows():
    temp = row['Compound_Name']
    full = temp
    offset = 0 
    for i in range(0, len(temp)-1):
        if temp[i+1] == 'C' and temp[i] != '=':
            new_cmps.append('{}O{}'.format(temp[:i+1], temp[i+1:]))
            parents.append(index)
            full = '{}O{}'.format(full[:i+1+offset], full[i+1+offset:])
            offset += 1
            ids.append(13)
    if offset > 0 :
        new_cmps += [full]
        parents.append(index)
        ids.append(14)
        
to_add = {'Compound_Name': new_cmps,  'ID': ids, 'Parent': parents}
oxy = pd.DataFrame(to_add, columns = ['Compound_Name', 'ID', 'Parent'])
all_compounds = pd.concat([fun_addition, oxy], ignore_index = True)
all_compounds = all_compounds.drop(['Place', 'Loc'], axis=1)
vps = []
for i in tqdm(range(len(all_compounds))):
    vps.append(vapour_pressure(all_compounds.loc[i, 'Compound_Name']))
#for i in range(len(all_compounds)):
#    
#    all_compounds.loc[i, 'Vapour_Pressure'] = vapour_pressure(all_compounds.loc[i, 'Compound_Name'])
    if i % 1000 == 0:
        print(i)
        
to_add1 = {'Compound_Name':list(all_compounds.loc[:, 'Compound_Name']), 'Vapour_Pressure': vps}
to_merge1 = pd.DataFrame(to_add1, columns = ['Compound_Name', 'Vapour_Pressure'])
all_compounds = pd.merge(all_compounds, to_merge1)
pickle.dump(all_compounds, open('save_data/all_compounds.p', 'wb'))

