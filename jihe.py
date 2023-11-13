from Bio.PDB import *
import numpy as np
from scipy.optimize import line_search
from Bio.PDB.Polypeptide import *
from Bio.PDB import Atom

def load_pdb(pdb_file):
    atoms = []
    parser = PDBParser()
    structure = parser.get_structure('test', pdb_file)
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if isinstance(atom.parent, Residue.Residue):
                        atom_name = atom.name.strip()
                        res_name = residue.resname.strip()
                        chain_id = chain.id.strip()
                        res_num = residue.id[1]
                        x, y, z = atom.coord  # 将坐标存储为numpy数组
                        element = atom.element.strip().upper()
                        new_atom = Atom.Atom(atom_name, (x, y, z), element, res_name, res_num, "", "", chain_id)
                        atoms.append(new_atom)
    return atoms

# 1. 读取 pdb 文件
atoms = load_pdb("/root/1a00B00.pdb")

# 2. 构建距离矩阵
n_atoms = len(atoms)
dist_matrix = np.zeros((n_atoms, n_atoms))
for i in range(n_atoms):
    for j in range(i + 1, n_atoms):
        dist_matrix[i, j] = np.linalg.norm(np.array(atoms[i].get_coord()) - np.array(atoms[j].get_coord()))
        dist_matrix[j, i] = dist_matrix[i, j]

# 3. 进行原子坐标的优化
energy = 0

coords = np.array([a.get_coord() for a in atoms])
def energy_function(x):
    e = 0
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            d = np.linalg.norm(x[(i*3):((i+1)*3)] - x[(j*3):((j+1)*3)])
            e += (dist_matrix[i,j] - d)**2
    return 0.5 * e

def gradient(x):
    g = np.zeros(x.shape)
    for i in range(n_atoms):
        for j in range(i+1, n_atoms):
            d = np.linalg.norm(x[(i*3):((i+1)*3)] - x[(j*3):((j+1)*3)])
            g[(i*3):((i+1)*3)] += (d - dist_matrix[i,j]) * (x[(i*3):((i+1)*3)] - x[(j*3):((j+1)*3)])
            g[(j*3):((j+1)*3)] += (d - dist_matrix[i,j]) * (x[(j*3):((j+1)*3)] - x[(i*3):((i+1)*3)])
    return g

def conjugate_gradient(x0, f, dfdx, maxiter=500, tol=1e-6):
    x = x0.copy()
    fval = f(x)
    g = dfdx(x)
    direction = -g
    for i in range(maxiter):
        alpha = line_search(f, dfdx, x.ravel(), direction.ravel())[0]
        x += alpha * direction
        fold = fval
        fval = f(x)
        gold = g
        g = dfdx(x)
        if np.linalg.norm(gold) == 0:
            beta = 0
        else:
            beta = np.dot(g, g - gold) / np.dot(gold, gold)
        direction = -g + beta * direction
        if np.linalg.norm(g) < tol:
            break
    return x, fval

optimized_coords, energy = conjugate_gradient(coords.reshape(-1, ), energy_function, gradient)
optimized_coords = optimized_coords.reshape((-1, 3))

# 构建新的结构
residues = {}
for atom in atoms:
    if atom.parent is None:
        continue
    residue_id = (atom.parent.id, atom.get_full_id()[2])
    if residue_id not in residues:
        residue_name = three_to_one(atom.parent.resname)
        residues[residue_id] = Residue.Residue(residue_id[1], residue_name, residue_id[0])

    new_atom = Atom.Atom(atom.name, optimized_coords[atom.serial_number-1], atom.bfactor, atom.occupancy, atom.altloc, atom.fullname, atom.serial_number)

    residues[residue_id].add(new_atom)

structure = Structure.Structure("reconstructed")
model = Model.Model(0)
for residue_id in residues:
    model.add(residues[residue_id])
structure.add(model)

# 将重构结果保存为 pdb 文件
io = PDBIO()
io.set_structure(structure)
io.save("/root/jihe")
print("PDB文件重构完成，结果保存为reconstructed.pdb")