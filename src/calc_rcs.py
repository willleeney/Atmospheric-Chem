def change_rc(fun_g, c_id, chain): # function that adds on a funtional group to the rate constant chains
    num_c = 0
    if len(chain) == 1: # if its CH4 it doesnt commute with the rest
        chain += [fun_g]
        return chain
    for i in range(len(chain)):
        if chain[i][0] != '(' and chain[i][0] != ')': #counts the number of carbon atoms 
            num_c += 1
        if num_c == c_id: #when the place is found 
            chain[i][0] = 'S' # changes the primary to sec 
            chain[i] += [fun_g] # adds the functional group to that atom
            if num_c == 1: # if it was at the start then the next atom needs the primary to be replaced
                for j in range(len(chain[1])):
                    if chain[1][j] == 'P':
                        chain[1][j] = 'S'
                        break 
            else:
                for j in range(i, 0, -1): # if not at start then needs to backtrack to find the atom it was attached to
                    if chain[j][0] == chain[i][1]:
                        for k in range(len(chain[j])):
                            if chain[j][k] == 'P':
                                chain[j][k] = 'S'
                                break 

            return chain

def cal_rc(chain, f_const, f_group):
    r_const = 0 # rate constant values - can be changed 
    k_pri = 0.136
    k_sec=0.934
    k_ter=1.94
    fx=1.23
    fch3=1
    for i in range(len(chain)):
        mult = 0
        if chain[i][0] != '(' and chain[i][0] != ')': # on carbon atom
            mult = 1
            if chain[i][0] == 'P': #times by the K value
                mult *= k_pri
            elif chain[i][0] == 'S':
                mult *= k_sec
            elif chain[i][0] == 'T':
                mult *= k_ter
            for j in range(1, len(chain[i])): # then times by the rest of the F(X)
                if chain[i][j] == 'P':
                    mult *= fch3
                elif chain[i][j] == f_group:
                    mult *= f_const
                else: 
                    mult *= fx
        r_const += mult 

    return r_const

all_compounds = pickle.load(open('save_data/all_compounds.p', 'rb'))
start_compounds = pickle.load(open('save_data/start_compounds.p', 'rb'))

rate_cmps = start_compounds.iloc[0:150].copy() # og alkanes
rate_chain = []
for index, row in rate_cmps.iterrows():
    temp_cmp = row['Compound_Name']
    #print(temp_cmp)
    for i in range(len(temp_cmp)): # first identify carbon atom to type 
        if temp_cmp[i] == 'C':
            if i == 0:
                #print('P')
                rate_chain += [[['P']]]
            elif i == (len(temp_cmp) - 1):
                #print('P')
                rate_chain[index] += [['P']]
            elif temp_cmp[i+1] == ')':
                #print('P')
                rate_chain[index] += [['P']]
            elif temp_cmp[i+1] == 'C':
                #print('S')
                rate_chain[index] += [['S']]
            elif temp_cmp[i+1:i+5] == '(C)(' or temp_cmp[i+1:i+6] == '(CC)(':
                #print('Q')
                rate_chain[index] += [['Q']]
            else: 
                #print('T')
                rate_chain[index] += [['T']]
        else: 
            #print(temp_cmp[i])
            rate_chain[index] += [temp_cmp[i]] #else to keep track of brackets


for i in range(1, len(rate_chain)): # finds what each atom is attached to 
    cmp_chain = rate_chain[i]
    #print(cmp_chain)
    for j in range(len(cmp_chain)):
        brac_lvl = 0
        fn = 0
        inc = 0
        if j == 0: # if at start - then primary
            #print('1') # next cmp_chain 
            
            while(fn != 1): #fwd # next letter on level 
                if cmp_chain[j+1+inc] == '(':
                    brac_lvl += 1
                elif cmp_chain[j+1+inc] == ')':
                    brac_lvl -= 1
                elif brac_lvl == 0:
                    cmp_chain[j] += cmp_chain[j+1+inc]
                    fn = 1
                inc += 1
                
        elif cmp_chain[j] == ['P']:
            #print('2')
            
            if cmp_chain[j-1] != '(' and cmp_chain[j-1] != ')' : # if last is not bracket then just the on before
                cmp_chain[j] += cmp_chain[j-1][0] 
            else:
                if cmp_chain [j-1] == ')':  # if its next to a closed chain then find next on level
                    while(fn != 1): #bck # last letter on level 
                        if cmp_chain[j-1-inc] == '(':
                            brac_lvl -= 1
                        elif cmp_chain[j-1-inc] == ')':
                            brac_lvl += 1
                        elif brac_lvl == 0:
                            cmp_chain[j] += cmp_chain[j-1-inc][0]
                            fn = 1
                        inc += 1
                else: # if its in a chain by itself find the one on the next level down
                    while(fn != 1): #bck # last letter on level 
                        if cmp_chain[j-1-inc] == '(':
                            brac_lvl -= 1
                        elif cmp_chain[j-1-inc] == ')':
                            brac_lvl += 1
                        elif brac_lvl < 0:
                            cmp_chain[j] += cmp_chain[j-1-inc][0]
                            fn = 1
                        inc += 1
                    
                
        elif cmp_chain[j] == ['S']:
            #print('3')
            # finds the one before on its level 
            while(fn != 1): #bck # last letter on level 
                if cmp_chain[j-1-inc] == '(':
                    brac_lvl -= 1
                elif cmp_chain[j-1-inc] == ')':
                    brac_lvl += 1
                elif brac_lvl <= 0:
                    cmp_chain[j] += cmp_chain[j-1-inc][0]
                    fn = 1
                inc += 1
                
            cmp_chain[j] += cmp_chain[j+1] # finds the next one which will always be an atom because of defintion 
            
        elif cmp_chain[j] == ['T']:
            #print('4')
            
            while(fn != 1): #next
                if cmp_chain[j+1+inc] != '(' and cmp_chain[j+1+inc] != ')':
                    cmp_chain[j] += cmp_chain[j+1+inc]
                    fn = 1
                inc += 1
                
            fn = 0 
            inc = 0
            while(fn != 1): #bck # last letter on level 
                if cmp_chain[j-1-inc] == '(':
                    brac_lvl -= 1
                elif cmp_chain[j-1-inc] == ')':
                    brac_lvl += 1
                elif brac_lvl <= 0:
                    cmp_chain[j] += cmp_chain[j-1-inc][0]
                    fn = 1
                inc += 1
                
            fn = 0 
            inc = 0 
            while(fn != 1): #fwd # next letter on level 
                if cmp_chain[j+1+inc] == '(':
                    brac_lvl += 1
                elif cmp_chain[j+1+inc] == ')':
                    brac_lvl -= 1
                elif brac_lvl == 0:
                    cmp_chain[j] += cmp_chain[j+1+inc]
                    fn = 1
                inc += 1
            
        elif cmp_chain[j] == ['Q']:
            #print('5')
            
            while(fn != 1): #bck # last letter on level 
                if cmp_chain[j-1-inc] == '(':
                    brac_lvl -= 1
                elif cmp_chain[j-1-inc] == ')':
                    brac_lvl += 1
                elif brac_lvl <= 0:
                    cmp_chain[j] += cmp_chain[j-1-inc][0]
                    fn = 1
                inc += 1 

            fn = 0 
            inc = 0 
            while(fn != 1): #fwd # next letter on level 
                if cmp_chain[j+1+inc] == '(':
                    brac_lvl += 1
                elif cmp_chain[j+1+inc] == ')':
                    brac_lvl -= 1
                elif brac_lvl == 0:
                    cmp_chain[j] += cmp_chain[j+1+inc]
                    fn = 1
                inc += 1

            fn = 0 
            inc = 1 
            op = 0 
            while(fn != 1): #fwd # next letter on level 
                if cmp_chain[j+1+inc] == '(':
                    brac_lvl += 1
                elif cmp_chain[j+1+inc] == ')':
                    brac_lvl -= 1
                elif brac_lvl == 0:
                    cmp_chain[j] += cmp_chain[j+1+inc]
                    fn = 1
                inc += 1

            fn = 0
            inc = 0 
            while(fn != 1): #next
                if cmp_chain[j+1+inc] != '(' and cmp_chain[j+1+inc] != ')':
                    cmp_chain[j] += cmp_chain[j+1+inc]
                    fn = 1
                inc += 1

    
og_rates = [] # calculates og rates 
for i in range(len(rate_chain)):
    og_rates.append(cal_rc(rate_chain[i], 'W', 0)) # the W and 0 are aribitay options 

    
all_compounds_rate = all_compounds.copy()
all_compounds_rate['Rate_Constant'] = pd.Series(og_rates, index=rate_cmps.index)

    
rate_cmps2 = start_compounds.iloc[150:1175].copy() # forget about the double bonds rates - yep same
new_chain = []
new_cmps = []
i = 0
for index,row in rate_cmps2.iterrows(): # removes the double bond
    temp_cmp = row['Compound_Name']
    temp = []
    for j in range(len(temp_cmp)):
        if temp_cmp[j] != '=':
            temp += temp_cmp[j]
    new_cmps.append(''.join(temp))
    new_cmps[i] = Chem.MolToSmiles(Chem.MolFromSmiles(new_cmps[i]))
    i += 1
    
j = 150
for i in range(len(new_cmps)): # gets the rate constant for its parent 
    p = start_compounds[start_compounds['Compound_Name'] == new_cmps[i]].index[0]
    all_compounds_rate.at[j, 'Rate_Constant'] = all_compounds_rate.at[int(p), 'Rate_Constant']
    j += 1
# its parent isnt stored because for functional groups alkanes and alkenes are the parents 
    
fun_id = [1, 3, 4, 5, 7, 8, 9]
fun_const = [3.4, 8.8, 0.74, 0.14, 0.38, 0.3, 0.53]
f_start = ['Original','O', 'OO', 'C(=O)', 'OC(=O)', 'N#C', 'N', 'Cl', 'Br', 'I', 'OCOO', '(=O)(OCOO)', 'OOCOO', 'O(mid)', 'O(all)']
fun_key = ['Original','O', 'OO', 'C=O', 'C(=O)O', 'C#N', 'N', 'Cl', 'Br', 'I', 'OCOO', '(=O)(OCOO)', 'OOCOO', 'O(mid)', 'O(all)']

for lp in range(len(fun_id)):
    f_id = fun_id[lp] # gets id 
    f_group = fun_key[f_id] # gets fun. group
    f_len = len(f_group)
    rate_const = all_compounds[all_compounds['ID'] == f_id].copy() # get compounds with functional group id
    rate_exit = []
    rate_id = []
    for index, row in rate_const.iterrows():
        temp_cmp = row['Compound_Name']
        if temp_cmp[0:len(f_start[f_id])] == f_start[f_id]: # if fun. group is at the start some groups are different
            rate_exit.append(temp_cmp[len(f_start[f_id]):])
            rate_id.append(1)
        else:
            r_id = 0 # counts the number of Cs after the fun_groups appears
            for i in range(len(temp_cmp)):
                if temp_cmp[i:i+2+f_len] == '({})'.format(f_group):
                    rate_exit.append(temp_cmp[0:i]+temp_cmp[i+2+f_len:])
                    rate_id.append(r_id)
                elif temp_cmp[i] == 'C':
                    r_id += 1
                    
    new_chain = []
    for i in range(len(rate_exit)):
        temp_cmp = rate_exit[i]
        temp = []
        for j in range(len(temp_cmp)):
            if temp_cmp[j] != '=':
                temp += temp_cmp[j]
        rate_exit[i] = ''.join(temp)
        rate_exit[i] = Chem.MolToSmiles(Chem.MolFromSmiles(rate_exit[i]))
        the_id = all_compounds[all_compounds['Compound_Name'] == rate_exit[i]].index[0]
        new_chain.append(change_rc(f_group, rate_id[i], copy.deepcopy(rate_chain[the_id])))
        
    rates = [] #calculate the br rates
    index = rate_const.index
    for i in range(len(new_chain)):
        all_compounds_rate.at[index[i], 'Rate_Constant'] = cal_rc(new_chain[i], f_group, fun_const[lp])  

pickle.dump(all_compounds_rate, open('all_compounds_rate.p', 'wb'))