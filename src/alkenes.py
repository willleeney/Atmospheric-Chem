################################## FINDING ALLENES AND ALKENES #################################################
def new_loc(n, nest, brac):
    # INPUT: n-number of hydrocarbons, nest-whether to use nested brackets or not, brac-the unnested brackets
    # OUTPUT: nest-whether to use nested brackets or not, brac-the next unnested brackets
    
    # gives the new location 
    for k in range(0, len(brac)):
        if brac[k] != (n-(1+k)):
            brac[k] += 1
        else:
            if len(brac) == n - 1: 
                nest = 1
            elif k == len(brac) - 1:
                brac = np.append(brac, 1)
                
    return nest, brac


from rdkit import Chem

alkanes = pickle.load(open('alkanes.p', 'rb'))

allenes = {'C': 1, 'CC': 2, 'CCC': 3, 'CCCC': 4, 'CCCCC': 5, 'CCCCCC': 6, 'CCCCCCC': 7, 'CCCCCCCC': 8, 'CCCCCCCCC': 9, 'CCCCCCCCCC': 10}
allenes = pd.DataFrame({'Compound_Name': list(allenes)}, columns = ['Compound_Name'])

for i, row in alkanes.iterrows(): 
    equ = [1]
    nest = 0
    if i != 0:
        while(nest == 0):
            iso = row['Compound_Name']
            strn = len(iso)
            nest, equ = new_loc(strn, nest, equ)
            
            for i in range(0, len(equ)):
                iso = iso[:equ[len(equ)-(i+1)]+i] + '=' + iso[equ[len(equ)-(i+1)]+i:]
                
            if IsoValid(allenes, iso) == True:
                allenes = allenes.append({'Compound_Name': Chem.MolToSmiles(Chem.MolFromSmiles(iso))}, ignore_index=True)
    
    
  
            
########################################### DELETING ALLENES ######################################################


alkenes = pd.DataFrame({'Compound_Name': list(allenes)[1:]}, columns = ['Compound_Name'])


for i, row in allenes.iterrows():
    mol = row['Compound_Name']
    dont = 0
    for k in range(0, len(mol)):
        if mol[k] == '=' and k <= len(mol) - 3:
            if mol[k+2] == '=':
                dont = 1
    if dont == 0:
        alkenes = alkenes.append({'Compound_Name': Chem.MolToSmiles(Chem.MolFromSmiles(mol))}, ignore_index=True)

start_compounds = pd.merge(alkanes, alkenes, how = 'outer')
pickle.dump(start_compounds, open("start_compounds.p", "wb"))