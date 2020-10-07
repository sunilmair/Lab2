from labutil.src.plugins.pwscf import *
from ase.spacegroup import crystal
from ase.io import write
import matplotlib.pyplot as plt

import numpy as np


def make_struc(alat):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    # set primitive_cell=False if you want to create a simple cubic unit cell with 8 atoms
    gecell = crystal('Ge', [(0, 0, 0)], spacegroup=227, cellpar=[alat, alat, alat, 90, 90, 90], primitive_cell=True)
    # check how your cell looks like
    # write('s.cif', gecell)
    structure = Struc(ase2struc(gecell))
    return structure

def make_struc_problem3(alat, z_displacement):
    """
    Creates the crystal structure using ASE.
    :param alat: Lattice parameter in angstrom
    :return: structure object converted from ase
    """
    # set primitive_cell=False if you want to create a simple cubic unit cell with 8 atoms
    gecell = crystal('Ge', [(0, 0, 0)], spacegroup=227, cellpar=[alat, alat, alat, 90, 90, 90], primitive_cell=True)
    print(gecell.positions)
    gecell.positions[0] = (gecell.positions[0][0], gecell.positions[0][1], gecell.position[0][2] + z_displacement*alat)
    print(gecell.positions)
    # check how your cell looks like
    # write('s.cif', gecell)
    structure = Struc(ase2struc(gecell))
    return structure


def compute_energy(alat, nk, ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potname = 'Ge.pz-bhs.UPF'
    pseudopath = os.environ['ESPRESSO_PSEUDO']
    potpath = os.path.join(pseudopath, potname)
    pseudopots = {'Ge': PseudoPotential(name=potname, path=potpath, ptype='uspp', element='Ge', functional='LDA')}
    struc = make_struc(alat=alat)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab2/Problem3/t1", str(nk)))
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'scf',
            'pseudo_dir': pseudopath,
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none',
        },
        'SYSTEM': {
            'ecutwfc': ecut,
             },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.5,
            'conv_thr': 1e-7,
        },
        'IONS': {},
        'CELL': {},
        })

    output_file = run_qe_pwscf(runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts)
    output = parse_qe_pwscf_output(outfile=output_file)
    return output


def compute_energy_problem3(alat, z_displacement, nk, ecut):
    """
    Make an input template and select potential and structure, and the path where to run
    """
    potname = 'Ge.pz-bhs.UPF'
    pseudopath = os.environ['ESPRESSO_PSEUDO']
    potpath = os.path.join(pseudopath, potname)
    pseudopots = {'Ge': PseudoPotential(name=potname, path=potpath, ptype='uspp', element='Ge', functional='LDA')}
    struc = make_struc_problem3(alat=alat, z_displacement=z_displacement)
    kpts = Kpoints(gridsize=[nk, nk, nk], option='automatic', offset=False)
    runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab2/Problem3/t1", str(nk)))
    input_params = PWscf_inparam({
        'CONTROL': {
            'calculation': 'scf',
            'pseudo_dir': pseudopath,
            'outdir': runpath.path,
            'tstress': True,
            'tprnfor': True,
            'disk_io': 'none',
        },
        'SYSTEM': {
            'ecutwfc': ecut,
             },
        'ELECTRONS': {
            'diagonalization': 'david',
            'mixing_beta': 0.5,
            'conv_thr': 1e-7,
        },
        'IONS': {},
        'CELL': {},
        })

    output_file = run_qe_pwscf(runpath=runpath, struc=struc,  pseudopots=pseudopots,
                               params=input_params, kpoints=kpts)
    output = parse_qe_pwscf_output(outfile=output_file)
    return output

def lattice_scan():
    nk = 9
    ecut = 30.0
    alat = 5.0
    output = compute_energy(alat=alat, ecut=ecut, nk=nk)
    print(output)
    energy = output['energy']
    #print(energy)

def problem1():
    #change runpathline in compute_energy to:
    #runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab2/Problem1/t1", str(ecut)))

    nk = 4
    ecut_list = np.arange(5.0, 85.0, 5.0)
    alat = 5.0
    output = [compute_energy(alat=alat, ecut=ecut, nk=nk) for ecut in ecut_list]

def problem2():
    # change runpathline in compute_energy to:
    # runpath = Dir(path=os.path.join(os.environ['WORKDIR'], "Lab2/Problem2/t1", str(nk)))

    nk_list = np.arange(1, 10, 1)
    ecut = 10.0
    alat = 5.0
    output = [compute_energy(alat=alat, ecut=ecut, nk=nk) for nk in nk_list]

def problem3():

    nk = 4
    ecut_list = np.arange(5.0, 85.0, 5.0)
    alat = 5.0
    z_displacement = 0.05
    output = [compute_energy(alat=alat, z_displacement=z_displacement, ecut=ecut, nk=nk) for ecut in ecut_list]


if __name__ == '__main__':
    # put here the function that you actually want to run
    make_struc_problem3(5.0, 0.05) #check if struc is good
    #lattice_scan() #check if 'force' is in 'output'
    #problem3()
