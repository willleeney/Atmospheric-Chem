def nonest(n, nest, brac):
    # INPUT: n-number of hydrocarbons, nest-whether to use nested brackets or not, brac-the unnested brackets
    # OUTPUT: iso-the next isomer attempt, nest-whether to use nested brackets or not, brac-the unnested brackets

    # gives the new bracket location 
    nest, brac = new_loc(n, nest, brac)
            
    # creates a non nested isomer
    iso = 'C'
    for l in range(1,n):
        added = 0
        for i in range(0,len(brac)):
            if brac[i] == l:
                iso += '(C)'
                added = 1
        if added == 0:
            iso += 'C'
                  
    return iso, nest, brac

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


def new_path(n, path):
    # INPUT: n-number of carbons, path-the current location of the brackets # 
    # OUTPUT: path-the new location of the brackets #
    s = np.shape(path)
    
    for k in range(s[0], 0, -1): # for every row starting from the bottom
        
        if path[k-1][1] == n: # if the right is full
            if path[k-1][0] != (n - 1): # if the left isnt 
                path[k-1][0] += 1 # increment left
                path[k-1][1] = path[k-1][0] + 1 
                break
            elif k != 1: 
                path[k-1][0] = 1
                path[k-1][1] = n - k 
            else: # if all are full
                path = add_reset(s, path)
        else: 
            path[k-1][1] += 1 # if right not full increment
            break
        
    return path  

def new_Isomer(n, path): 
    # INPUT: n-number of carbons, path-the location of the brackets # 
    # OUTPUT: isomer-the isomer string #
    isomer = 'C'
    s = np.shape(path)
    
    for l in range(1, n+1): # for all Carbon atoms # 
        for i in range(0,s[0]):# for the number of brackets 
            if path[i][1] == l:
                isomer += ')'
            if path[i][0] == l: 
                isomer += '('
        isomer += 'C'  
        
    return isomer

def add_reset(s, path):
    # INPUT: s-current path shape path-current locations of brackets 
    # OUTPUT: path- the reset of path
    for l in range(0, s[0]):
        path[l][0] = 1
        path[l][1] = 2+l
    path = np.vstack([path, [1, s[0]+2]])
    
    return path

def IsoValid(Isomers, potential):
    # INPUT: Isomers-list of isomers for current No. Carbons, potential-potential new isomer # 
    # OUTPUT: Boolean validity of isomer #
    if Chem.MolFromSmiles(potential) is None: # if the SMILE string is not valid
        return False
    if Isomers.get(Chem.MolToSmiles(Chem.MolFromSmiles(potential))) != None:
        return False # then it is not valid  
    return True



############################################# FINDING ALKANES #####################################################

alkanes = {'C': 1, 'CC': 2, 'CCC': 3, 'CCCC': 4, 'CCCCC': 5, 'CCCCCC': 6, 'CCCCCCC': 7, 'CCCCCCCC': 8, 'CCCCCCCCC': 9, 'CCCCCCCCCC': 10}  # the list containing all isomers
no_isomers = [1,1,1,2,3,5,9,18,35,75] # the number of isomers for each hydrocarbon

for n in range(0,10): # the current hydrocarbon
    
    nest = 0
    path = [[1,1]] # the starting position of path
    brac = np.array([1]) 
    count = 1
    
    while (count != no_isomers[n]): # until all possible isomers are found for hydrocarbon
        if nest == 0:
            # finds a new isomer without nesting brackets 
            iso_try, nest, brac = nonest(n+1, nest, brac)
            
        elif nest == 1:
            # find new bracket location # 
            path = new_path(n, path)

            # find possible isomer 
            iso_try = new_Isomer(n, path)
        
        # check isomer is valid # 
        valid = IsoValid(alkanes, iso_try) 
        
        # if valid then added # 
        if valid == True:
            alkanes[Chem.MolToSmiles(Chem.MolFromSmiles(iso_try))] = n+1
            count += 1

pickle.dump(alkanes, open("alkanes.p", "wb"))