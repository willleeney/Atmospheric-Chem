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

def vapour_pressure(compound_name):
    # INPUT: compound_name: the SMILE string of a compound #
    # OUTPUT: v_p: the vapour pressure at 300 d Kelvin #
    data = { "compounds": [compound_name], "temperatures": [300], "vp_method": "nannoolal", "bp_method": "nannoolal"}
    headers = {'Accept' : 'application/json'}
    response = requests.post("http://umansysprop.seaes.manchester.ac.uk/api/vapour_pressure", data = json.dumps(data), headers=headers)
    r = response.json()
    
    return r[0]['data'][0]['value']