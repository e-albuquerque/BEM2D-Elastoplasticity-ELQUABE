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


def compute_plate_analytical_solution(coords, ge, xnu, traction_y=1.0):
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


def compute_cylinder_analytical_solution(coords, E, nu, a=1.0, b=1.5, pressure=10.0):
    """
    Computes the analytical 2D plane-strain Lamé solution for a thick-walled
    cylinder of inner radius `a` and outer radius `b` under internal pressure `P`
    (matching `cylinder_analytical_results.py`):
        u_r(r)       = ((1 + nu) * P * a^2) / (E * (b^2 - a^2)) * ((1 - 2*nu)*r + b^2 / r)
        sigma_r(r)   = (a^2 * P) / (b^2 - a^2) * (1 - b^2 / r^2)
        sigma_t(r)   = (a^2 * P) / (b^2 - a^2) * (1 + b^2 / r^2)
        sigma_vm(r)  = sqrt(sigma_r^2 - sigma_r * sigma_t + sigma_t^2)
    """
    x = coords[:, 0]
    y = coords[:, 1]
    r = np.hypot(x, y)
    factor = (a ** 2 * pressure) / (b ** 2 - a ** 2)
    sigma_r = factor * (1.0 - (b ** 2) / (r ** 2))
    sigma_theta = factor * (1.0 + (b ** 2) / (r ** 2))
    u_r = ((1.0 + nu) / E) * factor * ((1.0 - 2.0 * nu) * r + (b ** 2) / r)
    sigma_vm = np.sqrt(sigma_r ** 2 - sigma_r * sigma_theta + sigma_theta ** 2)
    return {
        'r': r,
        'u_r': u_r,
        'sigma_r': sigma_r,
        'sigma_theta': sigma_theta,
        'sigma_vm': sigma_vm,
    }


def cartesian_to_polar_results(coords, disp, stress):
    """
    Transforms Cartesian displacements (ux, uy) and stresses (sx, txy, sy, svm)
    into polar coordinates (r, u_r, sigma_r, sigma_theta, sigma_vm).
    """
    x = coords[:, 0]
    y = coords[:, 1]
    r = np.hypot(x, y)
    theta = np.arctan2(y, x)
    c = np.cos(theta)
    s = np.sin(theta)

    u_r = disp[:, 0] * c + disp[:, 1] * s
    sx = stress[:, 0]
    txy = stress[:, 1]
    sy = stress[:, 2]
    svm = stress[:, 3]

    sigma_r = sx * (c ** 2) + sy * (s ** 2) + 2.0 * txy * s * c
    sigma_theta = sx * (s ** 2) + sy * (c ** 2) - 2.0 * txy * s * c

    return {
        'r': r,
        'u_r': u_r,
        'sigma_r': sigma_r,
        'sigma_theta': sigma_theta,
        'sigma_vm': svm,
    }


def plot_validation(inp, out, inp_params=None, save_path='result_plot.png'):
    """
    Selects the analytical solution corresponding to the active problem in
    `input_data.py`, prints global relative L2 errors, and plots the validation
    comparison in `result_plot.png`.
    """
    if inp_params is None:
        inp_params = input_data()

    problem_type = inp_params.get('problem_type', 'plate')
    bou_coords = inp['bou_coords']
    int_coords = inp['int_coords']
    ge = inp['ge']
    xnu = inp['xnu']

    dpbou = out['dpbou']
    dpint = out['dpint']
    stbou = out['stbou']
    stint = out['stint']

    if problem_type == 'cylinder':
        params = inp_params.get('analytical_params', {})
        a = params.get('a', 1.0)
        b = params.get('b', 1.5)
        pressure = params.get('pressure', 10.0)
        E = inp_params['E']
        nu = inp_params['nu']

        r_cont = np.linspace(a, b, 100)
        cont_coords = np.column_stack([r_cont, np.zeros_like(r_cont)])
        exact_cont = compute_cylinder_analytical_solution(cont_coords, E, nu, a=a, b=b, pressure=pressure)

        exact_bou = compute_cylinder_analytical_solution(bou_coords, E, nu, a=a, b=b, pressure=pressure)
        exact_int = compute_cylinder_analytical_solution(int_coords, E, nu, a=a, b=b, pressure=pressure)

        num_bou = cartesian_to_polar_results(bou_coords, dpbou, stbou)
        num_int = cartesian_to_polar_results(int_coords, dpint, stint)

        all_ur_num = np.concatenate([num_bou['u_r'], num_int['u_r']])
        all_ur_ex = np.concatenate([exact_bou['u_r'], exact_int['u_r']])
        err_ur = np.linalg.norm(all_ur_num - all_ur_ex) / np.linalg.norm(all_ur_ex)

        all_sr_num = np.concatenate([num_bou['sigma_r'], num_int['sigma_r']])
        all_sr_ex = np.concatenate([exact_bou['sigma_r'], exact_int['sigma_r']])
        err_sr = np.linalg.norm(all_sr_num - all_sr_ex) / np.linalg.norm(all_sr_ex)

        all_st_num = np.concatenate([num_bou['sigma_theta'], num_int['sigma_theta']])
        all_st_ex = np.concatenate([exact_bou['sigma_theta'], exact_int['sigma_theta']])
        err_st = np.linalg.norm(all_st_num - all_st_ex) / np.linalg.norm(all_st_ex)

        all_vm_num = np.concatenate([stbou[:, 3], stint[:, 3]])
        all_vm_ex = np.concatenate([exact_bou['sigma_vm'], exact_int['sigma_vm']])
        err_vm = np.linalg.norm(all_vm_num - all_vm_ex) / np.linalg.norm(all_vm_ex)

        print("=" * 78)
        print("ELQUABE 2D ELASTOSTATIC VALIDATION SUMMARY (CYLINDER UNDER INTERNAL PRESSURE)")
        print("=" * 78)
        print(f"Relative L2 Error in Radial Displacement u_r        : {err_ur:.4e}")
        print(f"Relative L2 Error in Radial Stress sigma_r          : {err_sr:.4e}")
        print(f"Relative L2 Error in Hoop Stress sigma_theta        : {err_st:.4e}")
        print(f"Relative L2 Error in von Mises Stress (All Nodes)   : {err_vm:.4e}")
        print("=" * 78)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Left subplot: Radial displacement u_r vs radius r
        ax1.plot(r_cont, exact_cont['u_r'], 'r-', linewidth=2, label='Analytical $u_r(r)$')
        ax1.plot(num_bou['r'], num_bou['u_r'], 'bo', markersize=5, label='BEM Boundary Nodes')
        ax1.plot(num_int['r'], num_int['u_r'], 'gs', markersize=5, label='BEM Internal Points')
        ax1.set_xlabel('Radius $r$ (m)')
        ax1.set_ylabel('Radial Displacement $u_r$ (m)')
        ax1.set_title('Radial Displacement $u_r$ vs. $r$')
        ax1.legend(loc='best')
        ax1.grid(True, linestyle='--', alpha=0.5)

        # Right subplot: Radial, hoop, and von Mises stresses vs radius r
        ax2.plot(r_cont, exact_cont['sigma_r'], 'r-', linewidth=2, label='Analytical $\\sigma_r(r)$')
        ax2.plot(r_cont, exact_cont['sigma_theta'], 'k--', linewidth=2, label='Analytical $\\sigma_\\theta(r)$')
        ax2.plot(r_cont, exact_cont['sigma_vm'], 'm-.', linewidth=1.8, label='Analytical $\\sigma_{\\mathrm{vM}}(r)$')
        ax2.plot(num_bou['r'], num_bou['sigma_r'], 'bo', markersize=5, label='BEM Boundary $\\sigma_r$')
        ax2.plot(num_int['r'], num_int['sigma_r'], 'c^', markersize=5, label='BEM Internal $\\sigma_r$')
        ax2.plot(num_bou['r'], num_bou['sigma_theta'], 'rs', markersize=4, fillstyle='none', label='BEM Boundary $\\sigma_\\theta$')
        ax2.plot(num_int['r'], num_int['sigma_theta'], 'gv', markersize=5, label='BEM Internal $\\sigma_\\theta$')
        ax2.plot(num_int['r'], num_int['sigma_vm'], 'm+', markersize=6, label='BEM Internal $\\sigma_{\\mathrm{vM}}$')
        ax2.set_xlabel('Radius $r$ (m)')
        ax2.set_ylabel('Stress (Pa)')
        ax2.set_title('Radial, Hoop & von Mises Stresses vs. $r$')
        ax2.legend(loc='best', fontsize=8)
        ax2.grid(True, linestyle='--', alpha=0.5)

        fig.suptitle('ELQUABE Validation: Quarter Cylinder Under Internal Pressure', fontsize=13)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Validation plot saved to '{save_path}'.")
        plt.close(fig)

    else:
        params = inp_params.get('analytical_params', {})
        traction_y = params.get('traction_y', 1.0)

        y_cont = np.linspace(0.0, np.max(bou_coords[:, 1]), 100)
        cont_coords = np.column_stack([np.full_like(y_cont, 0.5), y_cont])
        exact_cont = compute_plate_analytical_solution(cont_coords, ge, xnu, traction_y=traction_y)

        exact_bou = compute_plate_analytical_solution(bou_coords, ge, xnu, traction_y=traction_y)
        exact_int = compute_plate_analytical_solution(int_coords, ge, xnu, traction_y=traction_y)

        # Compute relative L2 errors across all nodes (boundary + internal)
        all_uy_num = np.concatenate([dpbou[:, 1], dpint[:, 1]])
        all_uy_ex = np.concatenate([exact_bou['uy'], exact_int['uy']])
        err_uy = np.linalg.norm(all_uy_num - all_uy_ex) / np.linalg.norm(all_uy_ex)

        all_vm_num = np.concatenate([stbou[:, 3], stint[:, 3]])
        all_vm_ex = np.concatenate([exact_bou['sigma_vm'], exact_int['sigma_vm']])
        err_vm = np.linalg.norm(all_vm_num - all_vm_ex) / np.linalg.norm(all_vm_ex)

        print("=" * 78)
        print("ELQUABE 2D ELASTOSTATIC VALIDATION SUMMARY (PLATE UNDER UNIAXIAL TRACTION)")
        print("=" * 78)
        print(f"Relative L2 Error in Vertical Displacement u_y      : {err_uy:.4e}")
        print(f"Relative L2 Error in von Mises Stress (All Nodes)   : {err_vm:.4e}")
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
        ax2.plot(y_cont, exact_cont['sigma_vm'], 'r-', linewidth=2, label=f'Analytical $\\sigma_{{\\mathrm{{vM}}}} = \\sigma_y = {traction_y:.1f}$')
        ax2.plot(y_cont, exact_cont['sigma_x'], 'k--', linewidth=1.5, label='Analytical $\\sigma_x = \\tau_{xy} = 0.0$')
        ax2.plot(bou_coords[:, 1], stbou[:, 3], 'bo', markersize=7, label='BEM Boundary $\\sigma_{\\mathrm{vM}}$')
        ax2.plot(int_coords[:, 1], stint[:, 3], 'gs', markersize=7, label='BEM Internal $\\sigma_{\\mathrm{vM}}$')
        ax2.plot(int_coords[:, 1], stint[:, 0], 'c^', markersize=6, label='BEM Internal $\\sigma_x$')
        ax2.plot(int_coords[:, 1], stint[:, 1], 'mv', markersize=6, label='BEM Internal $\\tau_{xy}$')
        ax2.set_xlabel('Vertical Coordinate $y$ (m)')
        ax2.set_ylabel('Stress (Pa)')
        ax2.set_title('Boundary & Internal Stresses vs. $y$')
        ax2.set_ylim(-0.2, 1.2 * traction_y)
        ax2.legend(loc='best', fontsize=8.5)
        ax2.grid(True, linestyle='--', alpha=0.5)

        fig.suptitle('ELQUABE 2D Elastostatic Validation: Analytical vs. BEM', fontsize=13)
        plt.tight_layout()
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Validation plot saved to '{save_path}'.")
        plt.close(fig)


def export_vtk_and_contours(
    out,
    inp_params=None,
    vtk_path='resultado_elquabe.vtk',
    disp_png_path='displacement_map.png',
    vm_png_path='von_mises_map.png',
):
    """
    Maps the BEM boundary and internal displacements and stresses (including
    von Mises equivalent stress from column 4 of `stbou.dat` and `stint.dat`)
    onto the Gmsh quadratic triangle mesh, exports a VTK file, and plots the
    2D mesh contour maps for displacement and von Mises stress.
    """
    if inp_params is None:
        inp_params = input_data()
    computed_data = compute_inodes(inp_params['file_name'], inp_params['bound_cond'])

    coord = computed_data['coordinates']
    ordered_boundary_nodes = computed_data['ordered_boundary_nodes']
    tri6 = computed_data['triangles']
    int_nodes = computed_data['internal_nodes']

    n_total = coord.shape[0]
    disp_vec = np.zeros((n_total, 3), dtype=np.float64)
    umag_field = np.zeros(n_total, dtype=np.float64)
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
        umag_field[node_id] = dpbou[bou_idx, 2]
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
        umag_field[node_id] = dpint[k, 2]
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
            'Displacement_Magnitude': umag_field,
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

    # 2a. Plot Displacement Contour Map (|u| for cylinder, u_y for plate)
    problem_type = inp_params.get('problem_type', 'plate')
    if problem_type == 'cylinder':
        disp_plot_field = umag_field
        disp_label = 'Radial Displacement $u_r$ (m)'
        disp_title = 'Radial Displacement $u_r$'
    else:
        disp_plot_field = uy_field
        disp_label = 'Vertical Displacement $u_y$ (m)'
        disp_title = 'Vertical Displacement $u_y$'

    fig1, ax1 = plt.subplots(figsize=(6, 5))
    ax1.triplot(triang, color=(0.0, 0.0, 0.0), linewidth=0.4)
    cf1 = ax1.tricontourf(triang, disp_plot_field, levels=20, cmap=cm.jet)
    fig1.colorbar(cf1, ax=ax1, label=disp_label)
    ax1.set_title(disp_title)
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
    ax2.set_title('von Mises Equivalent Stress $\\sigma_{\\mathrm{vM}}$')
    ax2.set_xlabel('$x$ (m)')
    ax2.set_ylabel('$y$ (m)')
    ax2.set_aspect('equal')
    plt.savefig(vm_png_path, dpi=300, bbox_inches='tight')
    print(f"von Mises contour map '{vm_png_path}' saved successfully.")
    plt.close(fig2)


if __name__ == '__main__':
    inp_params = input_data()
    inp = read_elquabe_input('eq12cb.dat')
    out = read_elquabe_outputs()
    plot_validation(inp, out, inp_params=inp_params, save_path='result_plot.png')
    msh_name = inp_params['file_name']
    if os.path.exists(msh_name + '.msh') or os.path.exists(msh_name + '.geo'):
        export_vtk_and_contours(out, inp_params=inp_params)
