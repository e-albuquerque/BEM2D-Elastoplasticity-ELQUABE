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
    boundary node coordinates, and internal point coordinates.
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

    int_coords = np.zeros((n_int, 2), dtype=np.float64)
    if n_int > 0:
        int_vals = [float(v) for v in lines[3 + ne].split()]
        int_coords = np.array(int_vals, dtype=np.float64).reshape((n_int, 2))

    return {
        'title': title,
        'ne': ne,
        'n_bou': n_bou,
        'n_int': n_int,
        'ge': ge,
        'xnu': xnu,
        'bou_coords': bou_coords,
        'int_coords': int_coords,
    }


def read_elquabe_outputs():
    """
    Reads numerical results exported by ELQUABE:
      - `dpbou.dat`: Boundary displacements (ux, uy, |u|)
      - `dpint.dat`: Internal displacements (ux, uy, |u|)
      - `stbou.dat`: Boundary tractions (tx, ty)
      - `stint.dat`: Internal stresses (sigma_x, tau_xy, sigma_y)
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
        sigma_x = 0, tau_xy = 0, sigma_y = traction_y
    """
    e_plane_strain = 2.0 * ge / (1.0 - xnu)
    y = coords[:, 1]
    uy_exact = (traction_y / e_plane_strain) * y
    sigma_y_exact = np.full_like(y, traction_y)
    sigma_x_exact = np.zeros_like(y)
    tau_xy_exact = np.zeros_like(y)
    return {
        'uy': uy_exact,
        'sigma_x': sigma_x_exact,
        'tau_xy': tau_xy_exact,
        'sigma_y': sigma_y_exact,
    }


def plot_validation(inp, out, save_path='result_plot.png'):
    """
    Compares numerical BEM displacements and stresses against the analytical
    2D elastostatic solution and saves the plot to `result_plot.png`.
    """
    bou_coords = inp['bou_coords']
    int_coords = inp['int_coords']
    ge = inp['ge']
    xnu = inp['xnu']

    dpbou = out['dpbou']
    dpint = out['dpint']
    stint = out['stint']

    y_cont = np.linspace(0.0, np.max(bou_coords[:, 1]), 100)
    cont_coords = np.column_stack([np.full_like(y_cont, 0.5), y_cont])
    exact_cont = compute_analytical_solution(cont_coords, ge, xnu, traction_y=1.0)

    exact_bou = compute_analytical_solution(bou_coords, ge, xnu, traction_y=1.0)
    exact_int = compute_analytical_solution(int_coords, ge, xnu, traction_y=1.0)

    # Compute relative L2 errors
    all_uy_num = np.concatenate([dpbou[:, 1], dpint[:, 1]])
    all_uy_ex = np.concatenate([exact_bou['uy'], exact_int['uy']])
    err_uy = np.linalg.norm(all_uy_num - all_uy_ex) / np.linalg.norm(all_uy_ex)
    err_sy = np.linalg.norm(stint[:, 2] - exact_int['sigma_y']) / np.linalg.norm(exact_int['sigma_y'])

    print("=" * 72)
    print("ELQUABE 2D ELASTOSTATIC VALIDATION SUMMARY")
    print("=" * 72)
    print(f"Relative L2 Error in Vertical Displacement u_y : {err_uy:.4e}")
    print(f"Relative L2 Error in Vertical Stress sigma_y   : {err_sy:.4e}")
    print("-" * 72)
    print("Internal Points Comparison (y, u_y_num, u_y_exact, sigma_y_num, sigma_y_exact):")
    for k in range(inp['n_int']):
        print(
            f"  Point {k+1} ({int_coords[k,0]:.2f}, {int_coords[k,1]:.2f}): "
            f"u_y = {dpint[k,1]:.6f} (exact: {exact_int['uy'][k]:.6f}) | "
            f"sigma_y = {stint[k,2]:.6f} (exact: {exact_int['sigma_y'][k]:.6f})"
        )
    print("=" * 72)

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

    # Right subplot: Internal stresses vs y
    ax2.plot(y_cont, exact_cont['sigma_y'], 'r-', linewidth=2, label='Analytical $\\sigma_y = 1.0$')
    ax2.plot(y_cont, exact_cont['sigma_x'], 'k--', linewidth=1.5, label='Analytical $\\sigma_x = \\tau_{xy} = 0.0$')
    ax2.plot(int_coords[:, 1], stint[:, 2], 'bs', markersize=7, label='BEM $\\sigma_y$')
    ax2.plot(int_coords[:, 1], stint[:, 0], 'go', markersize=6, label='BEM $\\sigma_x$')
    ax2.plot(int_coords[:, 1], stint[:, 1], 'm^', markersize=6, label='BEM $\\tau_{xy}$')
    ax2.set_xlabel('Vertical Coordinate $y$ (m)')
    ax2.set_ylabel('Stress Components (Pa)')
    ax2.set_title('Internal Point Stresses vs. $y$')
    ax2.set_ylim(-0.2, 1.2)
    ax2.legend(loc='best')
    ax2.grid(True, linestyle='--', alpha=0.5)

    fig.suptitle('ELQUABE 2D Elastostatic Validation: Analytical vs. BEM', fontsize=13)
    plt.tight_layout()
    plt.savefig(save_path, dpi=300, bbox_inches='tight')
    print(f"Validation plot saved to '{save_path}'.")
    plt.close(fig)


def export_vtk_and_contour(out, vtk_path='resultado_elquabe.vtk', png_path='displacement_map.png'):
    """
    Maps the BEM boundary and internal displacements back onto the Gmsh quadratic
    triangle mesh (`placa.msh`), exports a VTK file, and plots the 2D contour map.
    """
    inp_params = input_data()
    computed_data = compute_inodes(inp_params['file_name'], inp_params['bound_cond'])

    coord = computed_data['coordinates']
    elements = computed_data['elements']
    tri6 = computed_data['triangles']
    int_nodes = computed_data['internal_nodes']

    n_total = coord.shape[0]
    disp_vec = np.zeros((n_total, 3), dtype=np.float64)
    uy_field = np.zeros(n_total, dtype=np.float64)

    # Map boundary node displacements in the same order written by preproc.py
    dpbou = out['dpbou']
    bou_idx = 0
    for elem in elements:
        for node_id in (elem[0], elem[2]):
            disp_vec[node_id, 0] = dpbou[bou_idx, 0]
            disp_vec[node_id, 1] = dpbou[bou_idx, 1]
            uy_field[node_id] = dpbou[bou_idx, 1]
            bou_idx += 1

    # Map internal node displacements
    dpint = out['dpint']
    for k, node_id in enumerate(int_nodes):
        disp_vec[node_id, 0] = dpint[k, 0]
        disp_vec[node_id, 1] = dpint[k, 1]
        uy_field[node_id] = dpint[k, 1]

    # 1. Export VTK with quadratic triangle6 elements
    cells = [('triangle6', tri6)]
    mesh_out = meshio.Mesh(
        coord,
        cells,
        point_data={
            'Displacement_Y': uy_field,
            'Displacement': disp_vec,
        },
    )
    meshio.write(vtk_path, mesh_out)
    print(f"VTK file '{vtk_path}' saved successfully.")

    # 2. Subdivide each 6-node triangle into 4 linear triangles for Matplotlib contour plot
    tri3_list = []
    for t in tri6:
        v0, v1, v2, v3, v4, v5 = t
        tri3_list.append([v0, v3, v5])
        tri3_list.append([v3, v1, v4])
        tri3_list.append([v5, v4, v2])
        tri3_list.append([v3, v4, v5])

    tri3_array = np.array(tri3_list)

    fig, ax = plt.subplots(figsize=(6, 5))
    triang = mtri.Triangulation(coord[:, 0], coord[:, 1], tri3_array)

    ax.triplot(triang, color=(0.0, 0.0, 0.0), linewidth=0.4)
    contourf = ax.tricontourf(triang, uy_field, levels=20, cmap=cm.jet)
    fig.colorbar(contourf, ax=ax, label='Vertical Displacement $u_y$ (m)')

    ax.set_title('Vertical Displacement $u_y$ Contour (Quadratic Elements)')
    ax.set_xlabel('$x$ (m)')
    ax.set_ylabel('$y$ (m)')
    ax.set_aspect('equal')

    plt.savefig(png_path, dpi=300, bbox_inches='tight')
    print(f"Contour map '{png_path}' saved successfully.")
    plt.close(fig)


if __name__ == '__main__':
    inp = read_elquabe_input('eq12cb.dat')
    out = read_elquabe_outputs()
    plot_validation(inp, out, save_path='result_plot.png')
    if os.path.exists('placa.msh') or os.path.exists('placa.geo'):
        export_vtk_and_contour(out)
