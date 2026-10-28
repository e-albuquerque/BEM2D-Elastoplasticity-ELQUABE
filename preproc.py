import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import meshio
from input_data import input_data


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
    1-based BEM boundary connectivity (`conec`), domain internal cells (`triangle6`
    with 6 nodes or `quad8` with 8 nodes), 1-based BEM cell connectivity (`conec_cells`),
    internal nodes, and boundary condition mappings.
    """
    ensure_mesh_exists(file_name)
    mesh = meshio.read(file_name + '.msh')
    coordinates = mesh.points
    line_elements = mesh.cells_dict['line3']
    segments = mesh.cell_data_dict['gmsh:physical']['line3']
    edges = mesh.field_data

    # Detect quadratic domain cells: 8-node quadrilaterals ('quad8') or 6-node triangles ('triangle6')
    if 'quad8' in mesh.cells_dict:
        domain_cells = mesh.cells_dict['quad8']
        nodc = 8
    elif 'triangle6' in mesh.cells_dict:
        domain_cells = mesh.cells_dict['triangle6']
        nodc = 6
    else:
        raise ValueError("Mesh must contain 'triangle6' (6-node) or 'quad8' (8-node) internal cells.")

    triangle_elements = mesh.cells_dict.get('triangle6', domain_cells)

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

    # Build ordered list of boundary nodes and 1-based BEM connectivity matrix
    # Note: Gmsh line3 format is [node_start, node_end, node_mid].
    # BEM local node numbering is [node_start (xi=-1), node_mid (xi=0), node_end (xi=+1)].
    ordered_boundary_nodes = []
    gmsh_to_bem = {}
    for elem in line_elements:
        for gmsh_node in (int(elem[0]), int(elem[2]), int(elem[1])):
            if gmsh_node not in gmsh_to_bem:
                ordered_boundary_nodes.append(gmsh_node)
                gmsh_to_bem[gmsh_node] = len(ordered_boundary_nodes)

    conec = np.zeros((line_elements.shape[0], 3), dtype=np.int64)
    for i, elem in enumerate(line_elements):
        conec[i, 0] = gmsh_to_bem[int(elem[0])]  # Node 1 (xi = -1)
        conec[i, 1] = gmsh_to_bem[int(elem[2])]  # Node 2 (xi =  0, mid-node)
        conec[i, 2] = gmsh_to_bem[int(elem[1])]  # Node 3 (xi = +1)

    boundary_nodes = np.unique(line_elements)
    all_nodes = np.unique(domain_cells)
    internal_nodes = [int(node) for node in all_nodes if node not in boundary_nodes]

    # Map internal nodes to 1-based BEM global indices N + 1 .. N + L
    n_bou = len(ordered_boundary_nodes)
    for k, gmsh_node in enumerate(internal_nodes):
        gmsh_to_bem[int(gmsh_node)] = n_bou + k + 1

    # Build 1-based BEM internal cell connectivity matrix (NCELL x NODC)
    # For triangle6: [corner1, corner2, corner3, mid12, mid23, mid31]
    # For quad8    : [corner1, corner2, corner3, corner4, mid12, mid23, mid34, mid41]
    ncell = domain_cells.shape[0]
    conec_cells = np.zeros((ncell, nodc), dtype=np.int64)
    for icell, cell in enumerate(domain_cells):
        for j in range(nodc):
            conec_cells[icell, j] = gmsh_to_bem[int(cell[j])]

    return {
        'nodes': boundary_nodes,
        'ordered_boundary_nodes': ordered_boundary_nodes,
        'conec': conec,
        'internal_nodes': internal_nodes,
        'all_nodes': all_nodes,
        'coordinates': coordinates,
        'elements': line_elements,
        'triangles': triangle_elements,
        'domain_cells': domain_cells,
        'nodc': nodc,
        'conec_cells': conec_cells,
        'gmsh_to_bem': gmsh_to_bem,
        'segments': segments,
        'bc_info': bc_info,
    }


def mount_bcs(segments, bc_info, elements=None, coordinates=None):
    """
    Constructs boundary condition type (`KODE`) and value (`DFI`) arrays
    for each 3-node quadratic boundary element (6 DOFs per element).

    Supports:
      * type_x, type_y in (0, 1): Prescribed displacement (0) or traction (1)
        in global (x, y) coordinates.
      * type_x == 2: Prescribed normal traction/pressure (value_x = t_n) along
        the boundary, which is decomposed at each quadratic element node into
        global (t_x, t_y) traction components (type = 1 in both x and y).
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

                if tp_x == 2:
                    if elements is None or coordinates is None:
                        raise ValueError(
                            "elements and coordinates are required in mount_bcs when type_x == 2."
                        )
                    # Local BEM node ordering: node1 (xi=-1), node2 (xi=0), node3 (xi=+1)
                    n1, n3, n2 = int(elements[i, 0]), int(elements[i, 1]), int(elements[i, 2])
                    x1, y1 = coordinates[n1, 0], coordinates[n1, 1]
                    x2, y2 = coordinates[n2, 0], coordinates[n2, 1]
                    x3, y3 = coordinates[n3, 0], coordinates[n3, 1]

                    for j, xi in enumerate((-1.0, 0.0, 1.0)):
                        dn1 = xi - 0.5
                        dn2 = -2.0 * xi
                        dn3 = xi + 0.5

                        dxdxi = dn1 * x1 + dn2 * x2 + dn3 * x3
                        dydxi = dn1 * y1 + dn2 * y2 + dn3 * y3
                        jac = np.hypot(dxdxi, dydxi)

                        # Outward unit normal vector (nx, ny)
                        nx = dydxi / jac
                        ny = -dxdxi / jac

                        # Decompose normal traction vl_x into global x and y tractions
                        type_bcs[i, 2 * j] = 1
                        type_bcs[i, 2 * j + 1] = 1
                        val_bcs[i, 2 * j] = vl_x * nx
                        val_bcs[i, 2 * j + 1] = vl_x * ny
                else:
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


def show_geometry(coord, xplot, yplot, elem, tri6, type_bcs, val_bcs, save_path='geometry.png'):
    """
    Visualizes the quadratic mesh and global (x, y) boundary conditions of the 2D
    elastoplastic problem, saving the figure to `geometry.png`.
    """
    num_elements = elem.shape[0]
    b_nodes = np.zeros((3 * num_elements, 2), dtype=np.float64)
    bcs = np.zeros((3 * num_elements, 4), dtype=np.float64)

    for i in range(num_elements):
        # Local nodes in BEM order: [node1 (xi=-1), node2 (xi=0), node3 (xi=+1)]
        local_gmsh_nodes = [elem[i, 0], elem[i, 2], elem[i, 1]]
        for j, gnode in enumerate(local_gmsh_nodes):
            idx = 3 * i + j
            b_nodes[idx, 0] = coord[gnode, 0]
            b_nodes[idx, 1] = coord[gnode, 1]
            bcs[idx, 0] = type_bcs[i, 2 * j]      # bc_type in x (0=displacement, 1=traction)
            bcs[idx, 1] = val_bcs[i, 2 * j]       # bc_value in x
            bcs[idx, 2] = type_bcs[i, 2 * j + 1]  # bc_type in y (0=displacement, 1=traction)
            bcs[idx, 3] = val_bcs[i, 2 * j + 1]   # bc_value in y

    # Extract indices of element nodes with known displacement (0) and traction (1)
    indknown_ux = np.where(bcs[:, 0] == 0)[0]  # Dirichlet BC in x
    indknown_tx = np.where(bcs[:, 0] == 1)[0]  # Neumann BC in x
    indknown_uy = np.where(bcs[:, 2] == 0)[0]  # Dirichlet BC in y
    indknown_ty = np.where(bcs[:, 2] == 1)[0]  # Neumann BC in y

    # Determine plot limits to provide a 15% margin around the geometry
    x_min, x_max = np.min(coord[:, 0]), np.max(coord[:, 0])
    y_min, y_max = np.min(coord[:, 1]), np.max(coord[:, 1])
    delta_x, delta_y = x_max - x_min, y_max - y_min
    max_delta = max(delta_x, delta_y) if max(delta_x, delta_y) > 0 else 1.0
    max_bc_val = np.max(np.abs(bcs[:, [1, 3]]))
    arrow_scale = (10.0 * max_bc_val / max_delta) if max_bc_val > 0 else 10.0

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.grid(True, linestyle='--', alpha=0.4)

    # Plot the underlying cell mesh
    if tri6.shape[1] == 6:
        tri3 = tri6[:, :3]
        ax.triplot(coord[:, 0], coord[:, 1], tri3, color="black", linewidth=0.4)
    elif tri6.shape[1] == 8:
        for q in tri6:
            pts_q = coord[[q[0], q[4], q[1], q[5], q[2], q[6], q[3], q[7], q[0]], :2]
            ax.plot(pts_q[:, 0], pts_q[:, 1], color="black", linewidth=0.4)

    # Plot quadratic boundary element curves
    for i in range(num_elements):
        line_segments = LineCollection(
            [list(zip(xplot[i, :], yplot[i, :]))],
            linestyles='solid',
            colors='black',
            linewidths=0.8,
        )
        ax.add_collection(line_segments)

    # Plot global x-direction vectors (scaled by BC value) for known traction (blue)
    if len(indknown_tx) > 0:
        ax.quiver(
            b_nodes[indknown_tx, 0], b_nodes[indknown_tx, 1],
            bcs[indknown_tx, 1], np.zeros_like(bcs[indknown_tx, 1]),
            color="blue", angles='xy', scale_units='xy', scale=arrow_scale,
            width=0.004, headaxislength=0,
        )

    # Plot global x-direction vectors (scaled by BC value) for known displacement (red)
    if len(indknown_ux) > 0:
        ax.quiver(
            b_nodes[indknown_ux, 0], b_nodes[indknown_ux, 1],
            bcs[indknown_ux, 1], np.zeros_like(bcs[indknown_ux, 1]),
            color="red", angles='xy', scale_units='xy', scale=arrow_scale,
            width=0.004, headaxislength=0,
        )

    # Plot global y-direction vectors (scaled by BC value) for known traction (blue)
    if len(indknown_ty) > 0:
        ax.quiver(
            b_nodes[indknown_ty, 0], b_nodes[indknown_ty, 1],
            np.zeros_like(bcs[indknown_ty, 3]), bcs[indknown_ty, 3],
            color="blue", angles='xy', scale_units='xy', scale=arrow_scale,
            width=0.004, headaxislength=0,
        )

    # Plot global y-direction vectors (scaled by BC value) for known displacement (red)
    if len(indknown_uy) > 0:
        ax.quiver(
            b_nodes[indknown_uy, 0], b_nodes[indknown_uy, 1],
            np.zeros_like(bcs[indknown_uy, 3]), bcs[indknown_uy, 3],
            color="red", angles='xy', scale_units='xy', scale=arrow_scale,
            width=0.004, headaxislength=0,
        )

    # Plot nodes with known traction (blue markers) and known displacement (red markers)
    ax.plot(b_nodes[indknown_tx, 0], b_nodes[indknown_tx, 1], 'bx', markersize=6, label='Known $t_x$')
    ax.plot(b_nodes[indknown_ty, 0], b_nodes[indknown_ty, 1], 'bP', markersize=4, label='Known $t_y$')
    ax.plot(b_nodes[indknown_ux, 0], b_nodes[indknown_ux, 1], 'rx', markersize=6, label='Known $u_x$')
    ax.plot(b_nodes[indknown_uy, 0], b_nodes[indknown_uy, 1], 'rP', markersize=4, label='Known $u_y$')

    ax.set_aspect("equal")
    ax.set_xlim(x_min - 0.15 * delta_x, x_max + 0.15 * delta_x)
    ax.set_ylim(y_min - 0.15 * delta_y, y_max + 0.15 * delta_y)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('Mesh, Internal Cells and Boundary Conditions (ELQUABE)')

    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Geometry plot saved to '{save_path}'.")
    plt.close(fig)


def format_coord(val):
    val_rounded = round(float(val), 8)
    if abs(val_rounded) < 1e-12:
        return "0."
    if val_rounded.is_integer():
        return f"{int(val_rounded)}."
    return f"{val_rounded:.8f}".rstrip('0')


def format_bc_val(val):
    val_rounded = round(float(val), 8)
    if abs(val_rounded) < 1e-12:
        return "0.0"
    if val_rounded.is_integer():
        return f"{val_rounded:.1f}"
    return f"{val_rounded:.8f}".rstrip('0')


def create_input_file(filename, computed_data, int_points, type_bcs, val_bcs, input_params):
    """
    Writes the formatted input file (`eq12cb.dat`) for the Fortran ELQUABE 2D
    Elastoplastic BEM solver, including:
      1. Title
      2. Basic parameters (NE, L, GE, XNU)
      3. Boundary node coordinates
      4. Boundary element connectivity (CONEC)
      5. Boundary conditions (KODE, DFI)
      6. Internal point coordinates
      7. Internal cell parameters (NCELL, NODC, NSIG, TOLGP, NAUTO)
      8. Internal cell connectivity (LNDC: NCELL lines with NODC nodes each)
      9. Plasticity material properties (NCRIT, UNIAX, FRICT, HARDI, HARDK, EM)
     10. Incremental-iterative control parameters (NINCS, FINCR, MITER, TOLER)
    """
    with open(filename, 'w') as f:
        # 1. Title Line
        f.write(f"{input_params['title']}\n")

        nodes = computed_data['coordinates']
        elements = computed_data['elements']
        ordered_boundary_nodes = computed_data['ordered_boundary_nodes']
        conec = computed_data['conec']
        conec_cells = computed_data['conec_cells']
        nodc = computed_data['nodc']
        ncell = conec_cells.shape[0]
        num_elements = elements.shape[0]
        num_internal_points = len(int_points)

        shear_modulus = input_params['E'] / (2.0 * (1.0 + input_params['nu']))
        poisson_ratio = input_params['nu']

        # 2. Basic Parameter Line: NE, L, GE, XNU
        f.write(f" {num_elements},{num_internal_points},        {shear_modulus:.12g},{poisson_ratio:.4f}\n")

        # 3. Boundary Nodes Coordinates Line
        node_coords_list = []
        for gmsh_node_id in ordered_boundary_nodes:
            x_val, y_val = nodes[gmsh_node_id][0], nodes[gmsh_node_id][1]
            node_coords_list.append(f"{format_coord(x_val)}  {format_coord(y_val)}")

        node_coords_str = "  ".join(node_coords_list)
        f.write("  " + node_coords_str + "\n")

        # 4. Element Connectivity Matrix (NE lines: NO1 NO2 NO3 for each quadratic element)
        for i in range(num_elements):
            f.write(f" {conec[i, 0]} {conec[i, 1]} {conec[i, 2]}\n")

        # 5. Boundary Conditions Lines (6 pairs of KODE, DFI per element)
        for i in range(num_elements):
            f.write(f"{type_bcs[i, 0]} {format_bc_val(val_bcs[i, 0])} ")
            f.write(f"{type_bcs[i, 1]} {format_bc_val(val_bcs[i, 1])} ")
            f.write(f"{type_bcs[i, 2]} {format_bc_val(val_bcs[i, 2])} ")
            f.write(f"{type_bcs[i, 3]} {format_bc_val(val_bcs[i, 3])} ")
            f.write(f"{type_bcs[i, 4]} {format_bc_val(val_bcs[i, 4])} ")
            f.write(f"{type_bcs[i, 5]} {format_bc_val(val_bcs[i, 5])} \n")

        # 6. Internal Points Coordinates Line
        if num_internal_points > 0:
            int_points_str = " ".join(
                f"{format_coord(nodes[p][0])}  {format_coord(nodes[p][1])}" for p in int_points
            )
            f.write(" " + int_points_str + "\n")

        # 7. Internal Cell Control Parameters: NCELL, NODC, NSIG, TOLGP, NAUTO
        plast = input_params.get('plasticity', {})
        nsig = plast.get('nsig', 4)
        tolgp = plast.get('tolgp', 1.0e-5)
        nauto = plast.get('nauto', -6)
        f.write(f" {ncell} {nodc} {nsig} {tolgp:.6e} {nauto}\n")

        # 8. Internal Cell Connectivity Matrix (NCELL lines: NODC global BEM node indices per cell)
        for icell in range(ncell):
            cell_nodes_str = " ".join(str(int(nd)) for nd in conec_cells[icell, :])
            f.write(f" {cell_nodes_str}\n")

        # 9. Material Plasticity Properties: NCRIT, UNIAX, FRICT, HARDI, HARDK, EM
        ncrit = plast.get('ncrit', 2)
        uniax = plast.get('uniax', 24.0)
        frict = plast.get('frict', 0.0)
        hardi = plast.get('hardi', 0.0)
        hardk = plast.get('hardk', 0.0)
        em = plast.get('em', 1.0)
        f.write(f" {ncrit} {uniax:.8e} {frict:.8e} {hardi:.8e} {hardk:.8e} {em:.8e}\n")

        # 10. Incremental-Iterative Parameters: NINCS, FINCR, MITER, TOLER
        nincs = plast.get('nincs', 0)
        fincr = plast.get('fincr', 0.95)
        miter = plast.get('miter', 50)
        toler = plast.get('toler', 1.0e-5)
        f.write(f" {nincs} {fincr:.6f} {miter} {toler:.6e}\n")

    print(
        f"Input file '{filename}' generated successfully "
        f"(NE={num_elements}, L={num_internal_points}, NCELL={ncell}, NODC={nodc})."
    )


if __name__ == '__main__':
    inp_data = input_data()
    computed_data = compute_inodes(inp_data['file_name'], inp_data['bound_cond'])
    type_bcs, val_bcs = mount_bcs(
        computed_data['segments'],
        computed_data['bc_info'],
        computed_data['elements'],
        computed_data['coordinates'],
    )
    xplot, yplot = comp_xyplot(computed_data['elements'], computed_data['coordinates'])
    show_geometry(
        computed_data['coordinates'],
        xplot,
        yplot,
        computed_data['elements'],
        computed_data['domain_cells'],
        type_bcs,
        val_bcs,
    )

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
