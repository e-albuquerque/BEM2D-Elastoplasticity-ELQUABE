import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import meshio


def input_data():
    """
    Defines input parameters for the 2D Elastostatic Boundary Element Method (BEM)
    simulation using quadratic boundary elements (ELQUABE).

    Returns:
        dict: Dictionary containing:
            * title (str): Problem description title
            * bound_cond (dict): Boundary conditions for each physical edge
              (type: 0 = prescribed displacement, 1 = prescribed traction)
            * E (float): Young's modulus
            * nu (float): Poisson's ratio
            * file_name (str): Base name of the Gmsh .msh file
    """
    bound_cond = {
        'botton': {'type_x': 0, 'value_x': 0.0, 'type_y': 0, 'value_y': 0.0},
        'right':  {'type_x': 1, 'value_x': 0.0, 'type_y': 1, 'value_y': 0.0},
        'top':    {'type_x': 1, 'value_x': 0.0, 'type_y': 1, 'value_y': 1.0},
        'left':   {'type_x': 1, 'value_x': 0.0, 'type_y': 1, 'value_y': 0.0},
    }
    E = 1.0
    nu = 0.0
    file_name = 'placa'
    title = 'SQUARE PLATE UNDER UNIAXIAL STATIC TRACTION (4 QUADRATIC ELEMENTS)'

    return {
        'title': title,
        'bound_cond': bound_cond,
        'E': E,
        'nu': nu,
        'file_name': file_name,
    }


def ensure_mesh_exists(file_name):
    """
    Ensures `<file_name>.msh` exists. If missing, attempts to run Gmsh on
    `<file_name>.geo`, or generates the reference 4-element quadratic mesh
    directly via meshio.
    """
    msh_file = file_name + '.msh'
    if os.path.exists(msh_file):
        return

    geo_file = file_name + '.geo'
    if os.path.exists(geo_file):
        try:
            subprocess.run(
                ['gmsh', '-2', geo_file, '-o', msh_file, '-order', '2'],
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
            if os.path.exists(msh_file):
                return
        except (FileNotFoundError, subprocess.CalledProcessError):
            pass

    # Fallback: write the standard 13-node quadratic mesh for the unit square plate
    points = np.array([
        [0.0,  0.0,  0.0],  # 0: bottom-left corner
        [1.0,  0.0,  0.0],  # 1: bottom-right corner
        [1.0,  1.0,  0.0],  # 2: top-right corner
        [0.0,  1.0,  0.0],  # 3: top-left corner
        [0.5,  0.0,  0.0],  # 4: bottom mid-node
        [1.0,  0.5,  0.0],  # 5: right mid-node
        [0.5,  1.0,  0.0],  # 6: top mid-node
        [0.0,  0.5,  0.0],  # 7: left mid-node
        [0.5,  0.5,  0.0],  # 8: center internal node
        [0.75, 0.25, 0.0],  # 9: internal node
        [0.25, 0.25, 0.0],  # 10: internal node
        [0.25, 0.75, 0.0],  # 11: internal node
        [0.75, 0.75, 0.0],  # 12: internal node
    ], dtype=np.float64)

    line3 = np.array([
        [0, 1, 4],
        [1, 2, 5],
        [2, 3, 6],
        [3, 0, 7],
    ], dtype=np.int32)

    triangle6 = np.array([
        [0, 1, 8, 4,  9, 10],
        [3, 0, 8, 7, 10, 11],
        [1, 2, 8, 5, 12,  9],
        [2, 3, 8, 6, 11, 12],
    ], dtype=np.int32)

    cells = [('line3', line3), ('triangle6', triangle6)]
    cell_data = {
        'gmsh:physical': [
            np.array([1, 2, 3, 4], dtype=np.int32),
            np.array([5, 5, 5, 5], dtype=np.int32),
        ],
        'gmsh:geometrical': [
            np.array([1, 2, 3, 4], dtype=np.int32),
            np.array([1, 1, 1, 1], dtype=np.int32),
        ],
    }
    field_data = {
        'botton': np.array([1, 1], dtype=np.int32),
        'right':  np.array([2, 1], dtype=np.int32),
        'top':    np.array([3, 1], dtype=np.int32),
        'left':   np.array([4, 1], dtype=np.int32),
        'surf1':  np.array([5, 2], dtype=np.int32),
    }

    mesh = meshio.Mesh(points=points, cells=cells, cell_data=cell_data, field_data=field_data)
    meshio.write(msh_file, mesh, file_format='gmsh22')
    print(f"Mesh file '{msh_file}' generated.")


def compute_inodes(file_name, bound_cond):
    """
    Reads a Gmsh `.msh` file and extracts nodes, quadratic boundary elements,
    domain triangles, internal nodes, and boundary condition mappings.
    """
    ensure_mesh_exists(file_name)
    mesh = meshio.read(file_name + '.msh')
    coordinates = mesh.points
    line_elements = mesh.cells_dict['line3']
    triangle_elements = mesh.cells_dict['triangle6']
    segments = mesh.cell_data_dict['gmsh:physical']['line3']
    edges = mesh.field_data

    bc_info = {}
    for key, value in bound_cond.items():
        segment = mesh.field_data[key][0]
        bc_info[key] = {
            'type_x': value['type_x'],
            'value_x': value['value_x'],
            'segment': segment,
            'type_y': value['type_y'],
            'value_y': value['value_y'],
            'edge': edges[key],
        }

    boundary_nodes = np.unique(line_elements)
    all_nodes = np.unique(triangle_elements)
    internal_nodes = [node for node in all_nodes if node not in boundary_nodes]

    return {
        'nodes': boundary_nodes,
        'internal_nodes': internal_nodes,
        'all_nodes': all_nodes,
        'coordinates': coordinates,
        'elements': line_elements,
        'triangles': triangle_elements,
        'segments': segments,
        'bc_info': bc_info,
    }


def mount_bcs(segments, bc_info):
    """
    Constructs boundary condition type (`KODE`) and value (`DFI`) arrays
    for each 3-node quadratic boundary element (6 DOFs per element).
    """
    num_elements = len(segments)
    type_bcs = np.zeros((num_elements, 6), dtype=np.int64)
    val_bcs = np.zeros((num_elements, 6), dtype=np.float64)

    for i, segment_index in enumerate(segments):
        for _, bc_data in bc_info.items():
            if bc_data['segment'] == segment_index:
                tp_x = bc_data['type_x']
                vl_x = bc_data['value_x']
                tp_y = bc_data['type_y']
                vl_y = bc_data['value_y']
                type_bcs[i, :] = np.array([tp_x, tp_y, tp_x, tp_y, tp_x, tp_y])
                val_bcs[i, :] = np.array([vl_x, vl_y, vl_x, vl_y, vl_x, vl_y])
                break

    return type_bcs, val_bcs


def comp_xyplot(elem, nodes):
    """
    Interpolates coordinates along each quadratic boundary element for plotting.
    """
    npp = 8
    eta = np.linspace(-1.0, 1.0, npp)
    num_elements = elem.shape[0]

    xplot = np.zeros((num_elements, npp))
    yplot = np.zeros((num_elements, npp))

    for i in range(num_elements):
        node1, node3, node2 = elem[i]
        x1, y1 = nodes[node1][0:2]
        x2, y2 = nodes[node2][0:2]
        x3, y3 = nodes[node3][0:2]

        for k in range(npp):
            N1 = 0.5 * eta[k] * (eta[k] - 1.0)
            N2 = 1.0 - eta[k] ** 2
            N3 = 0.5 * eta[k] * (eta[k] + 1.0)
            xplot[i, k] = N1 * x1 + N2 * x2 + N3 * x3
            yplot[i, k] = N1 * y1 + N2 * y2 + N3 * y3

    return xplot, yplot


def show_geometry(nodes, xplot, yplot, elem, save_path='geometry.png'):
    """
    Plots the boundary element geometry and nodes, saving to `geometry.png`.
    """
    xnode = nodes[:, 0]
    ynode = nodes[:, 1]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.plot(xnode[elem[:, 0]], ynode[elem[:, 0]], 'gd', markersize=8, label='Element Extreme Nodes')
    ax.plot(xnode, ynode, 'ro', markersize=4, label='All Mesh Nodes')
    ax.set_aspect('equal')

    n_el = elem.shape[0]
    for i in range(n_el):
        x = xplot[i, :]
        y = yplot[i, :]
        line_segments = LineCollection([list(zip(x, y))], linestyles='solid', colors='tab:blue')
        ax.add_collection(line_segments)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Quadratic BEM Mesh & Internal Points (ELQUABE)')
    ax.legend(loc='best')
    ax.grid(True, linestyle='--', alpha=0.4)

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Geometry plot saved to '{save_path}'.")
    plt.close(fig)


def format_coord(val):
    val_rounded = round(val, 4)
    if val_rounded.is_integer():
        return f"{int(val_rounded)}."
    return f"{val_rounded}"


def create_input_file(filename, computed_data, int_points, type_bcs, val_bcs, input_params):
    """
    Writes the formatted input file (`eq12cb.dat`) for the Fortran ELQUABE solver.
    """
    with open(filename, 'w') as f:
        # 1. Title Line
        f.write(f"{input_params['title']}\n")

        nodes = computed_data['coordinates']
        elements = computed_data['elements']
        num_elements = elements.shape[0]
        num_internal_points = len(int_points)

        shear_modulus = input_params['E'] / (2.0 * (1.0 + input_params['nu']))
        poisson_ratio = input_params['nu']

        # 2. Basic Parameter Line: NE, L, GE, XNU
        f.write(f" {num_elements},{num_internal_points},        {shear_modulus:.1f},{poisson_ratio:.2f}\n")

        # 3. Boundary Nodes Coordinates Line (counter-clockwise: node 1 and mid-node per element)
        node_coords_list = []
        for elem in elements:
            no1_idx = elem[0]
            no_meio_idx = elem[2]
            x1, y1 = nodes[no1_idx][0], nodes[no1_idx][1]
            xm, ym = nodes[no_meio_idx][0], nodes[no_meio_idx][1]
            node_coords_list.append(f"{format_coord(x1)}  {format_coord(y1)}")
            node_coords_list.append(f"{format_coord(xm)}  {format_coord(ym)}")

        node_coords_str = "  ".join(node_coords_list)
        f.write("  " + node_coords_str + "\n")

        # 4. Boundary Conditions Lines (6 pairs of KODE, DFI per element)
        for i in range(num_elements):
            f.write(f"{type_bcs[i, 0]} {val_bcs[i, 0]:.1f} ")
            f.write(f"{type_bcs[i, 1]} {val_bcs[i, 1]:.1f} ")
            f.write(f"{type_bcs[i, 2]} {val_bcs[i, 2]:.1f} ")
            f.write(f"{type_bcs[i, 3]} {val_bcs[i, 3]:.1f} ")
            f.write(f"{type_bcs[i, 4]} {val_bcs[i, 4]:.1f} ")
            f.write(f"{type_bcs[i, 5]} {val_bcs[i, 5]:.1f} \n")

        # 5. Internal Points Coordinates Line
        int_points_str = " ".join(f"{nodes[p][0]:.1f}  {nodes[p][1]:.1f}" for p in int_points)
        f.write(" " + int_points_str + "\n")

    print(f"Input file '{filename}' generated successfully.")


if __name__ == '__main__':
    inp_data = input_data()
    computed_data = compute_inodes(inp_data['file_name'], inp_data['bound_cond'])
    type_bcs, val_bcs = mount_bcs(computed_data['segments'], computed_data['bc_info'])
    xplot, yplot = comp_xyplot(computed_data['elements'], computed_data['coordinates'])
    show_geometry(computed_data['coordinates'], xplot, yplot, computed_data['elements'])

    inp_data_filename = 'eq12cb.dat'
    int_points = computed_data['internal_nodes']
    create_input_file(
        inp_data_filename,
        computed_data,
        int_points,
        type_bcs,
        val_bcs,
        inp_data,
    )
