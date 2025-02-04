import os
import urllib.request
import urllib.parse
from urllib.error import URLError, HTTPError
import json
import re
import sys
import functools #implemented caching bc of high network usage. Things like FASTA sequences or ligand names usually don't change very often.

defaultCache='./pdbCache/'
pdbDatabase='https://files.rcsb.org/download/'
validFileTypes=['cif','pdb']

def url_exists(url:str):
    try:
        request = urllib.request.Request(url, method='HEAD')
        with urllib.request.urlopen(request) as response:
            return response.status in range(200, 400)
    except HTTPError as e:
        return False
    except URLError as e:
        return False
    

def checkInit(cache=defaultCache):
    if not os.path.isdir(cache): 
        os.mkdir(cache)

@functools.lru_cache(maxsize=50)
def get(target:str):
    checkInit()

    if os.path.exists(target):
        return(target)
    
    target=target.split('\\')[-1].split('/')[-1]
    type=target.split('.')
    if len(type) == 1:
        if not os.path.exists(defaultCache+target+'.cif'):
            url=pdbDatabase+target+'.cif'
            if url_exists(url):
                urllib.request.urlretrieve(url,defaultCache+target+'.cif')
            else:
                raise(ValueError,f'Error retrieving \'{target}\' from PDB: {url} could not be reached')
        targetPath=defaultCache+target+'.cif'
    
    elif type[-1] in validFileTypes:
        if not os.path.exists(defaultCache+target):
            url=pdbDatabase+target
            if url_exists(url):
                urllib.request.urlretrieve(url,defaultCache+target)
            else:
                raise(ValueError,f'Error retrieving \'{target}\' from PDB: {url} could not be reached')
        targetPath=defaultCache+target
    else:
        raise(ValueError,f'ERROR GETTING COORDINATE FILE \'{target}\': INVALID FILE TYPE.\nVALID FILE TYPES: {validFileTypes}')
    
    if os.path.exists(targetPath):
        return targetPath
    
@functools.lru_cache(maxsize=50)
def getFasta(upCode:str) -> str:
    url=f'https://rest.uniprot.org/uniprotkb/{upCode}.fasta'
    if url_exists(url):
        fasta=urllib.request.urlopen(url).read().decode('utf-8')
        seq=''
        parse=False
        for line in fasta:
            if line=='\n':
                parse=True
                continue
            if line == '>' and parse:
                break
            if parse:
                seq+=line.strip()

        return seq
    else:
        return None

# Filter for names. 
# Excludes []{} to avoid getting IUPAC names, and needs letters to avoid other codes
def valid_name(name:str) -> bool:
    if re.search(r'[\[\]{},]', name):
        return False
    if re.search('[a-zA-Z]', name) and re.search(r'^[^-]*-?[^-]*$', name) and re.search(r'^(?!.*\d{5,}).*$', name):
        return True
    else:
        return False

# Retrieves common name of a ligand (if possible)
@functools.lru_cache(maxsize=50)
def ccdName(ccd_code:str) -> str:
    '''
    Finds the common name for a chemical in the CCD. Returns the CCD code if a name could not be found.
    '''
    rcsb_url=f"https://data.rcsb.org/rest/v1/core/chemcomp/{ccd_code.upper()}" # First checks RCSB's accession of the CCD. This often has the best common name, but it isn't always there (ie. GW-4064)
    try:
        with urllib.request.urlopen(rcsb_url) as rcsb_response:
            rcsb_data = json.loads(rcsb_response.read().decode())
            synonyms = rcsb_data['rcsb_chem_comp_synonyms']
            for i in synonyms:
                if valid_name(i['name']):
                    return i['name']
            pubchem_id=None
            if 'rcsb_chem_comp_related' in rcsb_data:
                other_resources=rcsb_data['rcsb_chem_comp_related']
                for i in other_resources:
                    if i['related_mapping_method'] == 'matching InChIKey in PubChem': # for falling back on pubchem. All ligands in the CCD should have a corresponding pubchem ID.
                        pubchem_id = i['resource_accession_code']
            else: # sometimes the pubchem code isn't baked into the string, so an extra ping to pubchem is neccesary to get the right entry
                inchi = rcsb_data['rcsb_chem_comp_descriptor']['in_ch_i']
                
                base = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchi/cids/JSON"
                inchi_encoded = urllib.parse.quote(inchi)
                inchi_url = f"{base}?inchi={inchi_encoded}"
                try:
                    with urllib.request.urlopen(inchi_url) as inchi_response:
                        inchi_data=json.loads(inchi_response.read().decode())

                        if 'IdentifierList' in inchi_data and 'CID' in inchi_data['IdentifierList']:
                            pubchem_id = inchi_data['IdentifierList']['CID'][0]
                except HTTPError as e:
                    print(f"Error: {e.reason} - Unable to retrieve CID.",file=sys.stderr)
            if not pubchem_id:
                print(f'Could not find a common name for "{ccd_code}"', file=sys.stderr)
                return ccd_code
            pubchem_url=f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{pubchem_id}/synonyms/JSON" # Checks pubchem next. This usually has the name if the CCD doesn't.
            try:
                with urllib.request.urlopen(pubchem_url) as pchem_response:
                    pchem_data=json.loads(pchem_response.read().decode())
                    synonyms = pchem_data['InformationList']['Information'][0]['Synonym']
                    for name in synonyms:
                        if valid_name(name):
                            if 'CHEMBL' in name:
                                chembl_url=f"https://www.ebi.ac.uk/chembl/api/data/molecule/{name}.json" # As a last resort, looks up the CHEMBL entry since pubchem usually has this listed after the common ligand names. This usually doesn't give anything, but there are a couple of cases where it does.
                                try:
                                    with urllib.request.urlopen(chembl_url) as chembl_response:
                                        chembl_data=json.loads(chembl_response.read().decode())
                                        synonyms=chembl_data['molecule_synonyms']
                                        for i in synonyms:
                                            if all([True if value not in i['syn_type'] else False for value in ['CODE','NUMBER']]) and valid_name(i['molecule_synonym']):
                                                return i['molecule_synonym']
                                        break
                                except HTTPError as e:
                                    print(f'Error: {e.reason} - could not retrieve ligand \'{name}\' from CHEMBL.', file=sys.stderr)
                            else:
                                return name

            except HTTPError as e:
                print(f'Error: {e.reason} - could not retrieve ligand \'{pubchem_id}\' from PubChem.', file=sys.stderr)

    except HTTPError as e:
        print(f"Error: {e.reason} - could not retrieve ligand '{ccd_code}' from RCSB.", file=sys.stderr)

    print(f'Could not find a common name for "{ccd_code}"', file=sys.stderr) # Defaults back to CCD code if nothing could be found
    return ccd_code

@functools.lru_cache(maxsize=50)
def iupacName(ccd_code:str) -> str:
    '''
    Uses RCSB's API to get the IUPAC name for a chemical in the CCD
    '''
    rcsb_url=f"https://data.rcsb.org/rest/v1/core/chemcomp/{ccd_code.upper()}"
    try:
        with urllib.request.urlopen(rcsb_url) as rcsb_response:
            rcsb_data=json.loads(rcsb_response.read().decode())
            identifiers=rcsb_data['pdbx_chem_comp_identifier']
            for id in identifiers:
                if id['type']=='SYSTEMATIC NAME':
                    return id['identifier']
    except HTTPError as e:
        print(f"Error: {e.reason} - could not retrieve ligand '{ccd_code}' from RCSB.", file=sys.stderr)

@functools.lru_cache(maxsize=50)
def smiles(ccd_code:str) -> str:
    rcsb_url=f"https://data.rcsb.org/rest/v1/core/chemcomp/{ccd_code.upper()}"
    try:
        with urllib.request.urlopen(rcsb_url) as rcsb_response:
            rcsb_data=json.loads(rcsb_response.read().decode())
            smiles=rcsb_data['rcsb_chem_comp_descriptor']['smiles']
            return smiles
    except HTTPError as e:
        print(f"Error: {e.reason} - could not retrieve ligand '{ccd_code}' from RCSB.", file=sys.stderr)

@functools.lru_cache(maxsize=50)
def pubchem_diagram(smiles:str, fileName:str, base_resolution=500, scale_factor=20):
    
    atom_count = len(re.findall(r'[A-Z][a-z]?', smiles))
    
    png_size = min(base_resolution+atom_count*scale_factor,5000)
    
    
    formatted_smiles=urllib.parse.quote(smiles)
    pchem_url=f'https://https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{formatted_smiles}/PNG?image_size={png_size}'
    try:
        response = urllib.request.urlopen(pchem_url)
        if response.status == 200:
            with open(fileName,'wb') as png:
                png.write(response.read())
        else:
            print('ERROR fetching image')
    except Exception as e:
        print(f'Err: {e}')


def csvToList(csvPath:str):
    cl=open(csvPath,'r')
    lout=[]
    for line in cl:
        line=line.strip('\n').split(',')
        if len(line) == 1:
            line = line[0]
        lout.append(line)
    cl.close()
    return lout