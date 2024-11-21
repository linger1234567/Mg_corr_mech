import numpy as np
from ase.io import read
from ase.neighborlist import NeighborList
import time

# 定义原子半径
atomic_radii = {
    'H': 0.47,
    'O': 0.84,
    'Mg': 1.67
}

# 预计算成键距离
def precompute_bonding_distances(atomic_radii):
    bonding_distances = {}
    atom_types = list(atomic_radii.keys())
    for i, atom1 in enumerate(atom_types):
        for atom2 in atom_types[i:]:
            bonding_distances[(atom1, atom2)] = bonding_distances[(atom2, atom1)] = (atomic_radii[atom1] + atomic_radii[atom2]) * 0.90
    return bonding_distances

# 判定成键状态
def determine_bonding_status(atoms, bonding_distances_lookup):
    species_count = {}
    n = len(atoms)
    nl = NeighborList([2.6] * n, bothways=True, self_interaction=False)
    nl.update(atoms)
    for i in range(n):
        atom = atoms[i]
        symbol_i = atom.symbol
        neighbors = nl.get_neighbors(i)[0]
        MgO_count = 0
        MgH_count = 0
        OH_count = 0
        H2_count = 0
        for j in neighbors:
            distance = atoms.get_distance(i, j, True)
            bonding_distance = bonding_distances_lookup[(symbol_i, atoms[j].symbol)]
            if distance < bonding_distance:
                symbol_j = atoms[j].symbol
                if atom.symbol == 'Mg' and symbol_j == 'O':
                    MgO_count += 1
                elif atom.symbol == 'Mg' and symbol_j == 'H':
                    MgH_count += 1
                elif atom.symbol == 'O' and symbol_j == 'H':
                    OH_count += 1
                elif atom.symbol == 'H' and symbol_j == 'H':
                    H2_count += 1
        if atom.symbol == 'Mg':
            if MgO_count > 0:
                species_count[f'MgO{MgO_count}'] = species_count.get(f'MgO{MgO_count}', 0) + 1
            if MgH_count > 0:
                species_count[f'MgH{MgH_count}'] = species_count.get(f'MgH{MgH_count}', 0) + 1
    return species_count

# 主函数
def main(filename, start_frame=0, end_frame=None, step=5):
    start_time = time.time()
    all_atoms = read(filename, index=':')
    bonding_distances_lookup = precompute_bonding_distances(atomic_radii)

    with open('bonding_results.txt', 'w') as f:
        for frame_count, atoms in enumerate(all_atoms, 1):
            if start_frame <= frame_count <= end_frame or end_frame is None:
                if frame_count % step == 0:
                    results = determine_bonding_status(atoms, bonding_distances_lookup)
                    f.write(f"Frame {frame_count}:\n")
                    for species, count in results.items():
                        f.write(f'{species}: {count}\n')
                    print(f"已处理帧 {frame_count}")

    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"总用时: {elapsed_time} 秒")

# 调用主函数
if __name__ == "__main__":
    input_filename = 'dump.xyz'  # 替换为实际文件名
    start_frame = 1  # 开始处理的帧
    end_frame = 2000  # 结束处理的帧
    step = 20  # 每n帧抽取一帧
    main(input_filename, start_frame=start_frame, end_frame=end_frame, step=step)
