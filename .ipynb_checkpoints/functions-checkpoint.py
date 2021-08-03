import json

def getMolname(mol):
    with open('molnames.txt') as f:
        names = [el.split() for el in f.readlines()]
    return names[mol-1][1]


def openJSON(infile):
    try:
        with open(infile) as f:
            return json.load(f)
    except FileNotFoundError:
        return None
    
    
def getMult(mol):
    with open('magnetic.txt') as f:
        lines = [line.split() for line in f.readlines() if line.strip() != '']
    for line in lines:
        try:
            if int(line[0]) == mol:
                return int(float(line[4])) + 1
        except IndexError:
            return 1


def parseOutput(mol, func):
    fname = f'calcs/{mol}_{func}_.json'
    
    try:
        output = openJSON(fname)
        world_prec = output['input']['scf_calculation']['scf_solver']['final_prec']
        orbital_thrs = output['input']['scf_calculation']['scf_solver']['orbital_thrs']
        energy_thrs = output['input']['scf_calculation']['scf_solver']['energy_thrs']
        version = output['output']['provenance']['version']
        ncores = output['output']['provenance']['total_cores']
        walltime = output['output']['scf_calculation']['scf_solver']['wall_time']  # seconds
        energy = output['output']['properties']['scf_energy']['E_tot']
        mult = output['input']['molecule']['multiplicity']
        ncycles = len(output['output']['scf_calculation']['scf_solver']['cycles'])
        success = output['output']['success']

        row = [mol, getMolname(mol), func, mult, world_prec, orbital_thrs,
              energy_thrs, version, ncores, walltime, ncycles, energy, success]
    except Exception:
        print('Bad termination', fname)
        row = mol, getMolname(mol), func, None, None, None, None, None, None, None, None, None, False
        
    return row


def getSCFCycleData(mol, func):
    fname = f'calcs/{mol}_{func}_.json'
    with open(fname) as f:
        output = json.load(f)
        
    data = [(cycle['energy_total'], abs(cycle['mo_residual']), abs(cycle['energy_update'])) for cycle in output['output']['scf_calculation']['scf_solver']['cycles']]
    return data


def getThresholds(mol, func):
    fname = f'calcs/{mol}_{func}_.json'
    
    with open(fname) as f:
        output = json.load(f)
        
    orb = output['input']['scf_calculation']['scf_solver']['orbital_thrs']
    en = output['input']['scf_calculation']['scf_solver']['energy_thrs']
    return orb, en