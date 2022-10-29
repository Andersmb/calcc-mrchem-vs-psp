import json
import os
import pandas as pd

def getMolname(idx, proj):
    with open(f'/Users/abr121/Documents/dev/proj_psp/{proj}_molnames.txt') as f:
        names = [el.split() for el in f.readlines() if el.strip() != '']
    for line in names:
        if idx == line[0]:
            return line[1]


def openJSON(infile):
    try:
        with open(infile) as f:
            return json.load(f)
    except FileNotFoundError:
        return None
    
    
def getMult(mol, proj):
    with open(f'{proj}_magnetic.txt') as f:
        lines = [line.split() for line in f.readlines() if line.strip() != '']
    for line in lines:
        try:
            if line[0] == mol:
                if proj == 'row12':
                    return int(float(line[4])) + 1
                elif proj == 'row3':
                    return int(float(line[-1])) + 1
        except IndexError:
            return 1


def parseOutput(f, proj):
    try:
        atom, i, func, ext = os.path.basename(f).split('_')
        mol = '_'.join([atom, i])
    except ValueError:
        mol, func, ext = os.path.basename(f).split('_')

    try:
        output = openJSON(f)
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
        restriction = True if mult == 1 else False

        row = [mol, getMolname(mol, proj), func.upper(), mult, restriction, world_prec, orbital_thrs,
              energy_thrs, version, ncores, walltime, ncycles, energy, success]
        
    except Exception:
        print('Bad termination', f)
        row = mol, getMolname(mol), func.upper(), None, None, None, None, None, None, None, None, None, None, False
        
    return row


def getSCFCycleData(f):
    with open(f) as file:
        output = json.load(file)
        
    data = [(cycle['energy_total'], abs(cycle['mo_residual']), abs(cycle['energy_update'])) for cycle in output['output']['scf_calculation']['scf_solver']['cycles']]
    return data


def getThresholds(f):
    with open(f) as file:
        output = json.load(file)
        
    orb = output['input']['scf_calculation']['scf_solver']['orbital_thrs']
    en = output['input']['scf_calculation']['scf_solver']['energy_thrs']
    return orb, en


def parse_molecule(mol):
    "Wow, this is one ugly function!"
    components = []
    for i, char in enumerate(mol):
        if char.isupper():
            for j, subchar in enumerate(mol[i+1:]):
                if subchar.isupper():
                    components.append([mol[i:i+j+1]])
                    break
            else:
                components.append([mol[i:]])
    
    components = [el for sub in components for el in sub]
    for i, comp in enumerate(components):
        if comp[-1].isnumeric():
            components[i] = [comp[:-1] for _ in range(int(comp[-1]))]
            
    for i, comp in enumerate(components):
        if not isinstance(comp, list):
            components[i] = [comp]
            
    return [el for sub in components for el in sub]


def compute_atomization_energy(mol=None, func=None, project=None):
    total = pd.read_csv(f'{project}_total_energies.csv')
    atomic = pd.read_csv(f'{project}_atomic_energies.csv')
    e = total.loc[(total.Molecule == mol) & (total.Functional == func)].iloc[0].Energy
    for atom in parse_molecule(mol):
        e -= atomic.loc[(atomic.Functional == func) & (atomic.Atom == atom)].Allel.iloc[0]
    return e