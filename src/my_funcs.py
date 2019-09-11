def IsoValid(Isomers, potential):
    # INPUT: Isomers-list of isomers for current No. Carbons, potential-potential new isomer # 
    # OUTPUT: Boolean validity of isomer #
    
    if Chem.MolFromSmiles(potential) is None: # if the SMILE string is not valid
        return False
    iso = Chem.MolToSmiles(Chem.MolFromSmiles(potential))
    if type(Isomers) == dict:
        if Isomers.get(iso) != None:
            return False # then it is not valid  
    elif type(Isomers) == pd.core.frame.DataFrame:
        if Isomers[Isomers['Compound_Name'] == iso].empty == False:
            return False 
    return True
