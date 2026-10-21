import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import cm
import meshio

from preproc import input_data, compute_inodes


def read_elquabe_input(filename='eq12cb.dat'):
    """
    Reads the ELQUABE input file (`eq12cb.dat`) and returns material properties,
    boundary node coordinates, element connectivity (`conec`), and internal
    point coordinates.
    """
    with open(filename, 'r') as f:
        lines = [line.strip() for line in f if line.strip()]

    title = lines[0]
    params = lines[1].replace(',', ' ').split()
    ne = int(params[0])
    n_int = int(params[1])
    ge = float(params[2])
    xnu = float(params[3])

    n_bou = 2 * ne
    bou_vals = [float(v) for v in lines[2].split()]
    bou_coords = np.array(bou_vals, dtype=np.float64).reshape((n_bou, 2))

    # Read element connectivity matrix (ne lines: NO1 NO2 NO3)
    conec = np.zeros((ne, 3), dtype=np.int64)
    for i in range(ne):
        conec[i, :] = [int(v) for v in lines[3 + i].split()[:3]]

    int_coords = np.zeros((n_int, 2), dtype=np.float64)
    if n_int > 0:
        int_vals = [float(v) for v in lines[3 + 2 * ne].split()]
        int_coords = np.array(int_vals, dtype=np.float64).reshape((n_int, 2))

    return {
        'title': title,
        'ne': ne,
        'n_bou': n_bou,
        'n_int': n_int,
        'ge': ge,
        'xnu': xnu,
        'bou_coords': bou_coords,
        'conec': conec,
        'int_coords': int_coords,
    }


def read_elquabe_outputs():
    """
    Reads numerical results exported by ELQUABE:
      - `dpbou.dat`: Boundary displacements (ux, uy, |u|)
      - `dpint.dat`: Internal displacements (ux, uy, |u|)
      - `stbou.dat`: Boundary stresses (sigma_x, tau_xy, sigma_y, sigma_vm)
      - `stint.dat`: Internal stresses (sigma_x, tau_xy, sigma_y, sigma_vm)
    """
    dpbou = np.loadtxt('dpbou.dat')
    dpint = np.loadtxt('dpint.dat')
    stbou = np.loadtxt('stbou.dat')
    stint = np.loadtxt('stint.dat')
    return {
        'dpbou': dpbou,
        'dpint': dpint,
        'stbou': stbou,
        'stint': stint,
    }


def compute_analytical_solution(coords, ge, xnu, traction_y=1.0):
    """
    Computes the analytical 2D plane-strain elastostatic solution for a plate
    fixed at y = 0 and subjected to uniform vertical traction sigma_y = traction_y
    at the top edge:
        E' = 2 * G / (1 - nu)
        u_y(y) = (traction_y / E') * y
        sigma_x = 0, tau_xy = 0, sigma_y = traction_y, sigma_vm = |traction_y|
    """
    e_plane_strain = 2.0 * ge / (1.0 - xnu)
    y = coords[:, 1]
    uy_exact = (traction_y / e_plane_strain) * y
    sigma_y_exact = np.full_like(y, traction_y)
    sigma_x_exact = np.zeros_like(y)
    tau_xy_exact = np.zeros_like(y)
    sigma_vm_exact = np.full_like(y, np.abs(traction_y))
    return {
        'uy': uy_exact,
        'sigma_x': sigma_x_exact,
        'tau_xy': tau_xy_exact,
        'sigma_y': sigma_y_exact,
        'sigma_vm': sigma_vm_exact,
    }


def plot_validation(inp, out, save_path='result_plot.png'):
    """
    Compares numerical BEM displacements and stresses (including boundary stresses
    and von Mises equivalent stress) against the analytical 2D elastostatic solution
    and saves the plot to `result_plot.png`.
    """
    bou_coords = inp['bou_coords']
    int_coords = inp['int_coords']
    ge = inp['ge']
    xnu = inp['xnu']

    dpbou = out['dpbou']
    dpint = out['dpint']
    stbou = out['stbou']
    stint = out['stint']

    y_cont = np.linspace(0.0, np.max(bou_coords[:, 1]), 100)
    cont_coords = np.column_stack([np.full_like(y_cont, 0.5), y_cont])
    exact_cont = compute_analytical_solution(cont_coords, ge, xnu, traction_y=1.0)

    exact_bou = compute_analytical_solution(bou_coords, ge, xnu, traction_y=1.0)
    exact_int = compute_analytical_solution(int_coords, ge, xnu, traction_y=1.0)

    # Compute relative L2 errors across all nodes (boundary + internal)
    all_uy_num = np.concatenate([dpbou[:, 1], dpint[:, 1]])
    all_uy_ex = np.concatenate([exact_bou['uy'], exact_int['uy']])
    err_uy = np.linalg.norm(all_uy_num - all_uy_ex) / np.linalg.norm(all_uy_ex)

    all_vm_num = np.concatenate([stbou[:, 3], stint[:, 3]])
    all_vm_ex = np.concatenate([exact_bou['sigma_vm'], exact_int['sigma_vm']])
    err_vm = np.linalg.norm(all_vm_num - all_vm_ex) / np.linalg.norm(all_vm_ex)

    print("=" * 78)
    print("ELQUABE 2D ELASTOSTATIC VALIDATION SUMMARY")
    print("=" * 78)
    print(f"Relative L2 Error in Vertical Displacement u_y      : {err_uy:.4e}")
    print(f"Relative L2 Error in von Mises Stress (All Nodes)   : {err_vm:.4e}")
    print("-" * 78)
    print("Boundary Nodes von Mises Stress (x, y, sigma_x, tau_xy, sigma_y, sigma_vM):")
    for i in range(inp['n_bou']):
        print(
            f"  Node {i+1} ({bou_coords[i,0]:.2f}, {bou_coords[i,1]:.2f}): "
            f"sx={stbou[i,0]: .2e}, txy={stbou[i,1]: .2e}, "
            f"sy={stbou[i,2]:.6f}, vM={stbou[i,3]:.6f}"
        )
    print("-" * 78)
    print("Internal Points von Mises Stress (x, y, sigma_x, tau_xy, sigma_y, sigma_vM):")
    for k in range(inp['n_int']):
        print(
            f"  Point {k+1} ({int_coords[k,0]:.2f}, {int_coords[k,1]:.2f}): "
            f"sx={stint[k,0]: .2e}, txy={stint[k,1]: .2e}, "
            f"sy={stint[k,2]:.6f}, vM={stint[k,3]:.6f}"
        )
    print("=" * 78)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    # Left subplot: Vertical displacement u_y vs y
    ax1.plot(y_cont, exact_cont['uy'], 'r-', linewidth=2, label='Analytical Solution $u_y(y)$')
    ax1.plot(bou_coords[:, 1], dpbou[:, 1], 'bo', markersize=7, label='BEM Boundary Nodes')
    ax1.plot(int_coords[:, 1], dpint[:, 1], 'gs', markersize=7, label='BEM Internal Points')
    ax1.set_xlabel('Vertical Coordinate $y$ (m)')
    ax1.set_ylabel('Vertical Displacement $u_y$ (m)')
    ax1.set_title('Vertical Displacement $u_y$ vs. $y$')
    ax1.legend(loc='best')
    ax1.grid(True, linestyle='--', alpha=0.5)

    # Right subplot: Boundary & internal stresses (including von Mises) vs y
    ax2.plot(y_cont, exact_cont['sigma_vm'], 'r-', linewidth=2, label='Analytical $\\sigma_{\\mathrm{vM}} = \\sigma_y = 1.0$')
    ax2.plot(y_cont, exact_cont['sigma_x'], 'k--', linewidth=1.5, label='Analytical $\\sigma_x = \\tau_{xy} = 0.0$')
    ax2.plot(bou_coords[:, 1], stbou[:, 3], 'bo', markersize=7, label='BEM Boundary $\\sigma_{\\mathrm{vM}}$')
    ax2.plot(int_coords[:, 1], stint[:, 3], 'gs', markersize=7, label='BEM Internal $\\sigma_{\\mathrm{vM}}$')
    ax2.plot(int_coords[:, 1], stint[:, 0], 'c^', markersize=6, label='BEM Internal $\\sigma_x$')
    ax2.plot(int_coords[:, 1], stint[:, 1], 'mv', markersize=6, label='BEM Internal $\\tau_{xy}$')
    ax2.set_xlabel('Vertical Coordinate $y$ (m)')
    ax2.set_ylabel('Stress (Pa)')
    ax2.set_title('Boundary & Internal Stresses vs. $y$')
    ax2.set_ylim(-0.2, 1.2)
    ax2.legend(loc='best', fontsize=8.5)
    ax2.grid(True, linestyle='--', alpha=0.5)

    fig.suptitle('ELQUABE 2D Elastostatic Validation: Analytical vs. BEM', fontsize=13)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Validation plot saved to '{save_path}'.")
    plt.close(fig)


def export_vtk_and_contours(
    out,
    vtk_path='resultado_elquabe.vtk',
    disp_png_path='displacement_map.png',
    vm_png_path='von_mises_map.png',
):
    """
    Maps the BEM boundary and internal displacements and stresses (including
    von Mises equivalent stress from column 4 of `stbou.dat` and `stint.dat`)
    onto the Gmsh quadratic triangle mesh (`placa.msh`), exports a VTK file,
    and plots the 2D mesh contour maps for displacement and von Mises stress.
    """
    inp_params = input_data()
    computed_data = compute_inodes(inp_params['file_name'], inp_params['bound_cond'])

    coord = computed_data['coordinates']
    ordered_boundary_nodes = computed_data['ordered_boundary_nodes']
    tri6 = computed_data['triangles']
    int_nodes = computed_data['internal_nodes']

    n_total = coord.shape[0]
    disp_vec = np.zeros((n_total, 3), dtype=np.float64)
    uy_field = np.zeros(n_total, dtype=np.float64)
    sx_field = np.zeros(n_total, dtype=np.float64)
    txy_field = np.zeros(n_total, dtype=np.float64)
    sy_field = np.zeros(n_total, dtype=np.float64)
    vm_field = np.zeros(n_total, dtype=np.float64)

    # Map boundary node displacements (dpbou.dat) and stresses (stbou.dat)
    # using the exact ordered_boundary_nodes list from compute_inodes
    dpbou = out['dpbou']
    stbou = out['stbou']
    for bou_idx, node_id in enumerate(ordered_boundary_nodes):
        disp_vec[node_id, 0] = dpbou[bou_idx, 0]
        disp_vec[node_id, 1] = dpbou[bou_idx, 1]
        uy_field[node_id] = dpbou[bou_idx, 1]

        sx_field[node_id] = stbou[bou_idx, 0]
        txy_field[node_id] = stbou[bou_idx, 1]
        sy_field[node_id] = stbou[bou_idx, 2]
        vm_field[node_id] = stbou[bou_idx, 3]

    # Map internal node displacements (dpint.dat) and stresses (stint.dat)
    dpint = out['dpint']
    stint = out['stint']
    for k, node_id in enumerate(int_nodes):
        disp_vec[node_id, 0] = dpint[k, 0]
        disp_vec[node_id, 1] = dpint[k, 1]
        uy_field[node_id] = dpint[k, 1]

        sx_field[node_id] = stint[k, 0]
        txy_field[node_id] = stint[k, 1]
        sy_field[node_id] = stint[k, 2]
        vm_field[node_id] = stint[k, 3]

    # 1. Export VTK with quadratic triangle6 elements
    cells = [('triangle6', tri6)]
    mesh_out = meshio.Mesh(
        coord,
        cells,
        point_data={
            'Displacement': disp_vec,
            'Displacement_Y': uy_field,
            'Sigma_X': sx_field,
            'Tau_XY': txy_field,
            'Sigma_Y': sy_field,
            'Von_Mises_Stress': vm_field,
        },
    )
    meshio.write(vtk_path, mesh_out)
    print(f"VTK file '{vtk_path}' saved successfully.")

    # 2. Subdivide each 6-node triangle into 4 linear triangles for Matplotlib contour plots
    tri3_list = []
    for t in tri6:
        v0, v1, v2, v3, v4, v5 = t
        tri3_list.append([v0, v3, v5])
        tri3_list.append([v3, v1, v4])
        tri3_list.append([v5, v4, v2])
        tri3_list.append([v3, v4, v5])

    tri3_array = np.array(tri3_list)
    triang = mtri.Triangulation(coord[:, 0], coord[:, 1], tri3_array)

    # 2a. Plot Vertical Displacement Contour Map
    fig1, ax1 = plt.subplots(figsize=(6, 5))
    ax1.triplot(triang, color=(0.0, 0.0, 0.0), linewidth=0.4)
    cf1 = ax1.tricontourf(triang, uy_field, levels=20, cmap=cm.jet)
    fig1.colorbar(cf1, ax=ax1, label='Vertical Displacement $u_y$ (m)')
    ax1.set_title('Vertical Displacement $u_y$ (Quadratic Elements)')
    ax1.set_xlabel('$x$ (m)')
    ax1.set_ylabel('$y$ (m)')
    ax1.set_aspect('equal')
    plt.savefig(disp_png_path, dpi=300, bbox_inches='tight')
    print(f"Displacement contour map '{disp_png_path}' saved successfully.")
    plt.close(fig1)

    # 2b. Plot von Mises Equivalent Stress Contour Map with Mesh Overlay
    fig2, ax2 = plt.subplots(figsize=(6, 5))
    cf2 = ax2.tricontourf(triang, vm_field, levels=20, cmap=cm.jet)
    ax2.triplot(triang, color=(0.0, 0.0, 0.0), linewidth=0.5)
    ax2.plot(coord[:, 0], coord[:, 1], 'ko', markersize=3.5)
    fig2.colorbar(cf2, ax=ax2, label='von Mises Stress $\\sigma_{\\mathrm{vM}}$ (Pa)')
    ax2.set_title('von Mises Equivalent Stress $\\sigma_{\\mathrm{vM}}$ (Quadratic Elements)')
    ax2.set_xlabel('$x$ (m)')
    ax2.set_ylabel('$y$ (m)')
    ax2.set_aspect('equal')
    plt.savefig(vm_png_path, dpi=300, bbox_inches='tight')
    print(f"von Mises contour map '{vm_png_path}' saved successfully.")
    plt.close(fig2)


if __name__ == '__main__':
    inp = read_elquabe_input('eq12cb.dat')
    out = read_elquabe_outputs()
    plot_validation(inp, out, save_path='result_plot.png')
    if os.path.exists('placa.msh') or os.path.exists('placa.geo'):
        export_vtk_and_contours(out)
