################################## ADDING PARENT COMPOUND ID ######################################################
functional_groups = ['O', 'OO', 'C=O', 'C(=O)O', '#N', 'N', 'Cl', 'Br', 'I'] # each of the initial fun groups
ind = []
new_cmps = []
#ryb = ['rgb(1,1,1)']
#ryb += cl.scales['12']['qual']['Set3']
 ######## adds parent for each standard function group addition ######
 
for i in tqdm(range(len(start_compounds))): # all framework compounds
    temp_cmp = start_compounds.at[i, 'Compound_Name']
    for j in range(len(temp_cmp)): # each position in compound
        if (j == 0):
            for k in range(len(functional_groups)): # each functional group
                if k == 4: ##########################
                    cmp = '{}{}'.format('N#', temp_cmp[:]) # N bond is different because is a triple bond 
                    ind.append(int(i))
                    new_cmps.append(cmp)
                else: 
                    cmp = '{}{}'.format(functional_groups[k], temp_cmp[:])
                    ind.append(int(i))
                    new_cmps.append(cmp)
        elif(j == len(temp_cmp)-1):
            for k in range(len(functional_groups)): # each functional group
                cmp = '{}({}){}'.format(temp_cmp[:j+1], functional_groups[k], temp_cmp[j+1:])
                ind.append(int(i))
                new_cmps.append(cmp)
        elif temp_cmp[j+1] == ')':
            for k in range(len(functional_groups)): # each functional group
                cmp = '{}({}){}'.format(temp_cmp[:j+1], functional_groups[k], temp_cmp[j+1:])
                ind.append(int(i))
                new_cmps.append(cmp)
                
                
to_add = {'Compound_Name': new_cmps, 'Parent': ind}
merge = pd.DataFrame(to_add, columns = ['Compound_Name', 'Parent'])
fun_addition = pd.merge(fun_addition, merge, how='outer')#, left_on = ['Colour'], right_on = ['Vapour_Pressure'])

########### ASSIGNS INITIAL PARENT ID ###########
for i in range(1175):
    fun_addition.at[i, 'Parent'] = int(i)

############ CHANGES COLOUR FOR NEW FUN GROUPS ##########
for index, row in fun_addition.iterrows():
    temp = row['Compound_Name']
    col = row['Colour']
    if col == ryb[1]: 
        if 'COO' in temp: 
            fun_addition.at[index, 'Colour'] = ryb[10]
    if col == ryb[2]: # this one doesnt work 
        if 'OOCOO' in temp: 
            fun_addition.at[index, 'Colour'] = ryb[11]
    if col == ryb[4]: 
        if '(OCOO)' in temp:
            fun_addition.at[index, 'Colour'] = ryb[12]

 ######### ADDS PARENT MOLECULE FOR LONGER FUN GROUPS ####### this could be improved
example1 = fun_addition[fun_addition['Colour'] == ryb[1]].copy()
example2 = fun_addition[fun_addition['Colour'] == ryb[2]].copy()
example3 = fun_addition[fun_addition['Colour'] == ryb[4]].copy()
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

my_copy = fun_addition.copy()
i = 44411
for index, row in example1.iterrows():
    my_copy.at[i, 'Parent'] = row['Parent']
    i += 1
for index, row in example2.iterrows():
    my_copy.at[i, 'Parent'] = row['Parent']
    i += 1
for index, row in example3.iterrows():
    my_copy.at[i, 'Parent'] = row['Parent']
    i += 1
#addition = pd.merge(fun_addition, example1, how = 'outer') 
#addition = pd.merge(addition, example2, how = 'outer') 
#addition = pd.merge(addition, example3, how = 'outer') 