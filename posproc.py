import os
import numpy as np
from scipy.optimize import brentq
from scipy.integrate import cumulative_trapezoid
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import cm
import meshio

from preproc import input_data, compute_inodes


def read_elquabe_input(filename='eq12cb.dat'):
    """
    Reads the ELQUABE input file (`eq12cb.dat`) and returns material properties,
    boundary node coordinates, element connectivity (`conec`), internal point
    coordinates, and 2D elastoplastic cell/material data if present.
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
    next_line_idx = 3 + 2 * ne
    if n_int > 0:
        int_vals = [float(v) for v in lines[next_line_idx].split()]
        int_coords = np.array(int_vals, dtype=np.float64).reshape((n_int, 2))
        next_line_idx += 1

    # Read optional 2D elastoplasticity section
    ncell = 0
    nodc = 0
    nsig = 4
    tolgp = 1.0e-5
    nauto = -6
    lndc = None
    plastic_props = None

    if next_line_idx < len(lines):
        cell_header = lines[next_line_idx].replace(',', ' ').split()
        ncell = int(cell_header[0])
        nodc = int(cell_header[1])
        nsig = int(cell_header[2])
        tolgp = float(cell_header[3])
        nauto = int(cell_header[4])
        next_line_idx += 1

        if ncell > 0:
            lndc = np.zeros((ncell, nodc), dtype=np.int64)
            for i in range(ncell):
                lndc[i, :] = [int(v) for v in lines[next_line_idx + i].split()[:nodc]]
            next_line_idx += ncell

            mat_vals = lines[next_line_idx].replace(',', ' ').split()
            inc_vals = lines[next_line_idx + 1].replace(',', ' ').split()
            plastic_props = {
                'ncrit': int(mat_vals[0]),
                'uniax': float(mat_vals[1]),
                'frict': float(mat_vals[2]),
                'hardi': float(mat_vals[3]),
                'hardk': float(mat_vals[4]),
                'em': float(mat_vals[5]),
                'nincs': int(inc_vals[0]),
                'fincr': float(inc_vals[1]),
                'miter': int(inc_vals[2]),
                'toler': float(inc_vals[3]),
            }

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
        'ncell': ncell,
        'nodc': nodc,
        'nsig': nsig,
        'lndc': lndc,
        'plastic_props': plastic_props,
    }


def read_elquabe_outputs():
    """
    Reads numerical results exported by ELQUABE:
      - `dpbou.dat`: Boundary displacements (ux, uy, |u|)
      - `dpint.dat`: Internal displacements (ux, uy, |u|)
      - `stbou.dat`: Boundary stresses & plastic strain (sigma_x, tau_xy, sigma_y, sigma_vm, epstn)
      - `stint.dat`: Internal stresses & plastic strain (sigma_x, tau_xy, sigma_y, sigma_vm, epstn)
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
    at the top edge.
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
    Computes the analytical 2D plane-strain Lamé (elastic) solution for a
    thick-walled cylinder of inner radius `a` and outer radius `b` under
    internal pressure `P`.
    """
    x = coords[:, 0]
    y = coords[:, 1]
    r = np.hypot(x, y)
    factor = (a ** 2 * pressure) / (b ** 2 - a ** 2)
    sigma_r = factor * (1.0 - (b ** 2) / (r ** 2))
    sigma_theta = factor * (1.0 + (b ** 2) / (r ** 2))
    sigma_z = nu * (sigma_r + sigma_theta)
    u_r = ((1.0 + nu) / E) * factor * ((1.0 - 2.0 * nu) * r + (b ** 2) / r)
    sigma_vm = np.sqrt(
        0.5 * ((sigma_r - sigma_theta) ** 2 + (sigma_theta - sigma_z) ** 2 + (sigma_z - sigma_r) ** 2)
    )
    return {
        'r': r,
        'u_r': u_r,
        'sigma_r': sigma_r,
        'sigma_theta': sigma_theta,
        'sigma_vm': sigma_vm,
        'epstn': np.zeros_like(r),
    }


def compute_cylinder_elastoplastic_solution(coords, E, nu, sigma_y, a=1.0, b=1.5, pressure=10.0):
    """
    Computes the analytical 2D plane-strain elastoplastic solution (Lamé–Hill)
    for a thick-walled cylinder (a <= r <= b) under internal pressure P with
    von Mises yield stress sigma_y:
      - At the elastic-plastic interface r = c, the elastic outer region
        (c <= r <= b) satisfies the exact plane-strain von Mises yield condition:
            p_c(c) = sigma_y * (b^2 - c^2) / sqrt(3 * b^4 + (1 - 2*nu)^2 * c^4)
      - In the plastic region (a <= r <= c), radial equilibrium with the plane-strain
        von Mises yield condition gives:
            sigma_r(r)     = -p_c - (2 / sqrt(3)) * sigma_y * ln(c / r)
            sigma_theta(r) = sigma_r(r) + (2 / sqrt(3)) * sigma_y
        Setting sigma_r(a) = -P determines the elastic-plastic radius c.
      - Radial displacement u_r(r) is given by Lamé's exact solution in c <= r <= b
        and integrated inward across a <= r <= c using plastic incompressibility
        (eps_r^p + eps_theta^p + eps_z^p = 0) and plane strain (eps_z = 0):
            (1 / r) * d(r * u_r)/dr = ((1 + nu) * (1 - 2*nu) / E) * (sigma_r + sigma_theta)
      - Equivalent plastic strain epstn(r) in a <= r <= c is obtained from the
        total deviatoric hoop-radial strain minus the elastic deviatoric part:
            epstn(r) = (2 / sqrt(3)) * ((u_r / r - du_r / dr) / 2 - (1 + nu) / E * (sigma_theta - sigma_r) / 2)
    """
    tau_y = sigma_y / np.sqrt(3.0)

    def pc_func(c_val):
        return sigma_y * (b ** 2 - c_val ** 2) / np.sqrt(3.0 * (b ** 4) + ((1.0 - 2.0 * nu) ** 2) * (c_val ** 4))

    # Check if purely elastic
    p_yield_init = pc_func(a)
    if pressure <= p_yield_init:
        res = compute_cylinder_analytical_solution(coords, E, nu, a=a, b=b, pressure=pressure)
        res['c_radius'] = a
        return res

    # Solve for elastic-plastic interface radius c in (a, b)
    def residual_c(c_val):
        return pc_func(c_val) + 2.0 * tau_y * np.log(c_val / a) - pressure

    c_radius = brentq(residual_c, a, b)
    p_c = pc_func(c_radius)

    x = coords[:, 0]
    y = coords[:, 1]
    r = np.hypot(x, y)

    sigma_r = np.zeros_like(r)
    sigma_theta = np.zeros_like(r)
    sigma_vm = np.zeros_like(r)
    u_r = np.zeros_like(r)
    epstn = np.zeros_like(r)

    mask_pl = r <= c_radius
    mask_el = ~mask_pl

    # 1. Elastic region (c < r <= b)
    if np.any(mask_el):
        r_el = r[mask_el]
        factor_el = (c_radius ** 2 * p_c) / (b ** 2 - c_radius ** 2)
        sr_el = factor_el * (1.0 - (b ** 2) / (r_el ** 2))
        st_el = factor_el * (1.0 + (b ** 2) / (r_el ** 2))
        sz_el = nu * (sr_el + st_el)
        sigma_r[mask_el] = sr_el
        sigma_theta[mask_el] = st_el
        sigma_vm[mask_el] = np.sqrt(
            0.5 * ((sr_el - st_el) ** 2 + (st_el - sz_el) ** 2 + (sz_el - sr_el) ** 2)
        )
        u_r[mask_el] = ((1.0 + nu) / E) * factor_el * ((1.0 - 2.0 * nu) * r_el + (b ** 2) / r_el)
        epstn[mask_el] = 0.0

    # 2. Plastic region (a <= r <= c)
    if np.any(mask_pl):
        r_pl = r[mask_pl]
        sr_pl = -p_c - 2.0 * tau_y * np.log(c_radius / r_pl)
        st_pl = sr_pl + 2.0 * tau_y
        sigma_r[mask_pl] = sr_pl
        sigma_theta[mask_pl] = st_pl
        sigma_vm[mask_pl] = sigma_y

        # Radial displacement at r = c from elastic outer region
        factor_c = (c_radius ** 2 * p_c) / (b ** 2 - c_radius ** 2)
        u_c = ((1.0 + nu) / E) * factor_c * ((1.0 - 2.0 * nu) * c_radius + (b ** 2) / c_radius)

        # Exact closed-form integral of r * (sigma_r + sigma_theta) from r to c:
        # sigma_r + sigma_theta = 2 * (tau_y - p_c) - 4 * tau_y * ln(c / r)
        # int_r^c rho * (sigma_r(rho) + sigma_theta(rho)) d(rho)
        #   = (tau_y - p_c) * (c^2 - r^2) - 4 * tau_y * (0.25 * (c^2 - r^2) - 0.5 * r^2 * ln(c / r))
        #   = -p_c * (c^2 - r^2) + 2 * tau_y * r^2 * ln(c / r)
        k_vol = (1.0 + nu) * (1.0 - 2.0 * nu) / E
        integral_r_sum = -p_c * (c_radius ** 2 - r_pl ** 2) + 2.0 * tau_y * (r_pl ** 2) * np.log(c_radius / r_pl)
        r_ur_pl = c_radius * u_c - k_vol * integral_r_sum
        ur_pl = r_ur_pl / r_pl
        u_r[mask_pl] = ur_pl

        # Total radial strain du_r/dr = (1/r)*d(r*u_r)/dr - u_r/r
        dur_dr_pl = k_vol * (sr_pl + st_pl) - ur_pl / r_pl
        gamma_tot = (ur_pl / r_pl) - dur_dr_pl
        gamma_el = ((1.0 + nu) / E) * (st_pl - sr_pl)
        epstn[mask_pl] = np.maximum(0.0, (gamma_tot - gamma_el) / np.sqrt(3.0))

    return {
        'r': r,
        'u_r': u_r,
        'sigma_r': sigma_r,
        'sigma_theta': sigma_theta,
        'sigma_vm': sigma_vm,
        'epstn': epstn,
        'c_radius': c_radius,
        'p_c': p_c,
    }


def cartesian_to_polar_results(coords, disp, stress):
    """
    Transforms Cartesian displacements (ux, uy) and stresses
    (sx, txy, sy, svm[, epstn]) into polar coordinates
    (r, u_r, sigma_r, sigma_theta, sigma_vm, epstn).
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
    epstn = stress[:, 4] if stress.shape[1] >= 5 else np.zeros_like(svm)

    sigma_r = sx * (c ** 2) + sy * (s ** 2) + 2.0 * txy * s * c
    sigma_theta = sx * (s ** 2) + sy * (c ** 2) - 2.0 * txy * s * c

    return {
        'r': r,
        'u_r': u_r,
        'sigma_r': sigma_r,
        'sigma_theta': sigma_theta,
        'sigma_vm': svm,
        'epstn': epstn,
    }


def plot_validation(inp, out, inp_params=None, save_path='result_plot.png'):
    """
    Compares BEM numerical results against the analytical solution (elastoplastic
    Lamé–Hill solution when `ncell > 0`, or elastic Lamé solution when `ncell == 0`),
    prints global relative L2 errors, and saves the validation plot to `result_plot.png`.
    """
    if inp_params is None:
        inp_params = input_data()

    problem_type = inp_params.get('problem_type', 'cylinder')
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

        plast_cfg = inp_params.get('plasticity', {})
        sigma_y = plast_cfg.get('uniax', 24.0)
        if inp.get('plastic_props') is not None:
            sigma_y = inp['plastic_props']['uniax']

        r_cont = np.linspace(a, b, 200)
        cont_coords = np.column_stack([r_cont, np.zeros_like(r_cont)])

        # Purely elastic Lamé reference curve
        elastic_cont = compute_cylinder_analytical_solution(
            cont_coords, E, nu, a=a, b=b, pressure=pressure
        )

        # Elastoplastic analytical solution (Lamé–Hill)
        exact_cont = compute_cylinder_elastoplastic_solution(
            cont_coords, E, nu, sigma_y, a=a, b=b, pressure=pressure
        )
        exact_bou = compute_cylinder_elastoplastic_solution(
            bou_coords, E, nu, sigma_y, a=a, b=b, pressure=pressure
        )
        exact_int = compute_cylinder_elastoplastic_solution(
            int_coords, E, nu, sigma_y, a=a, b=b, pressure=pressure
        )
        c_radius = exact_cont['c_radius']

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

        all_vm_num = np.concatenate([num_bou['sigma_vm'], num_int['sigma_vm']])
        all_vm_ex = np.concatenate([exact_bou['sigma_vm'], exact_int['sigma_vm']])
        err_vm = np.linalg.norm(all_vm_num - all_vm_ex) / np.linalg.norm(all_vm_ex)

        all_ep_num = np.concatenate([num_bou['epstn'], num_int['epstn']])
        all_ep_ex = np.concatenate([exact_bou['epstn'], exact_int['epstn']])

        print("=" * 80)
        print("ELQUABE 2D ELASTOPLASTIC VALIDATION SUMMARY (CYLINDER UNDER INTERNAL PRESSURE)")
        print("=" * 80)
        print(f"Inner Radius a = {a:.3f} m, Outer Radius b = {b:.3f} m, Pressure P = {pressure:.2f} Pa")
        print(f"Yield Stress sigma_Y = {sigma_y:.2f} Pa, Analytical Plastic Radius c = {c_radius:.4f} m")
        print("-" * 80)
        print(f"Relative L2 Error in Radial Displacement u_r        : {err_ur:.4e} ({100.0*err_ur:.2f}%)")
        print(f"Relative L2 Error in Radial Stress sigma_r          : {err_sr:.4e} ({100.0*err_sr:.2f}%)")
        print(f"Relative L2 Error in Hoop Stress sigma_theta        : {err_st:.4e} ({100.0*err_st:.2f}%)")
        print(f"Relative L2 Error in von Mises Stress sigma_vM      : {err_vm:.4e} ({100.0*err_vm:.2f}%)")
        print(f"Max Equivalent Plastic Strain (BEM vs. Analytical)  : {np.max(all_ep_num):.4e} vs. {np.max(all_ep_ex):.4e}")
        print("=" * 80)

        fig, axes = plt.subplots(2, 2, figsize=(13, 10))
        ax1, ax2, ax3, ax4 = axes.flatten()

        # (a) Radial displacement u_r vs radius r
        ax1.plot(r_cont, elastic_cont['u_r'], 'k:', linewidth=1.5, label='Elastic Lamé $u_r(r)$')
        ax1.plot(r_cont, exact_cont['u_r'], 'r-', linewidth=2.2, label='Analytical Elastoplastic $u_r(r)$')
        ax1.plot(num_bou['r'], num_bou['u_r'], 'bo', markersize=4.5, alpha=0.8, label='BEM Boundary Nodes')
        ax1.plot(num_int['r'], num_int['u_r'], 'gs', markersize=4.0, alpha=0.7, label='BEM Internal Nodes')
        ax1.axvline(c_radius, color='darkorange', linestyle='--', linewidth=1.4, label=f'Plastic Radius $c = {c_radius:.3f}$ m')
        ax1.set_xlabel('Radius $r$ (m)')
        ax1.set_ylabel('Radial Displacement $u_r$ (m)')
        ax1.set_title('(a) Radial Displacement $u_r(r)$')
        ax1.legend(loc='best', fontsize=8.5)
        ax1.grid(True, linestyle='--', alpha=0.5)

        # (b) Radial & Hoop stresses vs radius r
        ax2.plot(r_cont, elastic_cont['sigma_theta'], 'k:', linewidth=1.4, label='Elastic Lamé $\\sigma_\\theta(r)$')
        ax2.plot(r_cont, exact_cont['sigma_r'], 'b-', linewidth=2.2, label='Analytical Elastoplastic $\\sigma_r(r)$')
        ax2.plot(r_cont, exact_cont['sigma_theta'], 'r-', linewidth=2.2, label='Analytical Elastoplastic $\\sigma_\\theta(r)$')
        ax2.plot(num_bou['r'], num_bou['sigma_r'], 'co', markersize=4.5, alpha=0.85, label='BEM Boundary $\\sigma_r$')
        ax2.plot(num_int['r'], num_int['sigma_r'], 'b^', markersize=3.8, alpha=0.65, label='BEM Internal $\\sigma_r$')
        ax2.plot(num_bou['r'], num_bou['sigma_theta'], 'ms', markersize=4.5, alpha=0.85, label='BEM Boundary $\\sigma_\\theta$')
        ax2.plot(num_int['r'], num_int['sigma_theta'], 'gv', markersize=3.8, alpha=0.65, label='BEM Internal $\\sigma_\\theta$')
        ax2.axvline(c_radius, color='darkorange', linestyle='--', linewidth=1.4, label=f'$c = {c_radius:.3f}$ m')
        ax2.set_xlabel('Radius $r$ (m)')
        ax2.set_ylabel('Stress (Pa)')
        ax2.set_title('(b) Radial Stress $\\sigma_r$ & Hoop Stress $\\sigma_\\theta$')
        ax2.legend(loc='best', fontsize=8)
        ax2.grid(True, linestyle='--', alpha=0.5)

        # (c) von Mises Equivalent Stress vs radius r
        ax3.plot(r_cont, elastic_cont['sigma_vm'], 'k:', linewidth=1.5, label='Elastic Lamé $\\sigma_{\\mathrm{vM}}(r)$')
        ax3.plot(r_cont, exact_cont['sigma_vm'], 'r-', linewidth=2.2, label='Analytical Elastoplastic $\\sigma_{\\mathrm{vM}}(r)$')
        ax3.plot(num_bou['r'], num_bou['sigma_vm'], 'bo', markersize=4.5, alpha=0.8, label='BEM Boundary $\\sigma_{\\mathrm{vM}}$')
        ax3.plot(num_int['r'], num_int['sigma_vm'], 'gs', markersize=4.0, alpha=0.7, label='BEM Internal $\\sigma_{\\mathrm{vM}}$')
        ax3.axvline(c_radius, color='darkorange', linestyle='--', linewidth=1.4, label=f'Plastic Radius $c = {c_radius:.3f}$ m')
        ax3.set_xlabel('Radius $r$ (m)')
        ax3.set_ylabel('von Mises Stress $\\sigma_{\\mathrm{vM}}$ (Pa)')
        ax3.set_title('(c) Equivalent von Mises Stress $\\sigma_{\\mathrm{vM}}(r)$')
        ax3.legend(loc='best', fontsize=8.5)
        ax3.grid(True, linestyle='--', alpha=0.5)

        # (d) Equivalent Plastic Strain vs radius r
        ax4.plot(r_cont, exact_cont['epstn'], 'r-', linewidth=2.2, label='Analytical $\\bar{\\varepsilon}^p(r)$')
        ax4.plot(num_bou['r'], num_bou['epstn'], 'bo', markersize=4.5, alpha=0.8, label='BEM Boundary $\\bar{\\varepsilon}^p$')
        ax4.plot(num_int['r'], num_int['epstn'], 'gs', markersize=4.0, alpha=0.7, label='BEM Internal $\\bar{\\varepsilon}^p$')
        ax4.axvline(c_radius, color='darkorange', linestyle='--', linewidth=1.4, label=f'Plastic Radius $c = {c_radius:.3f}$ m')
        ax4.set_xlabel('Radius $r$ (m)')
        ax4.set_ylabel('Equivalent Plastic Strain $\\bar{\\varepsilon}^p$')
        ax4.set_title('(d) Equivalent Plastic Strain $\\bar{\\varepsilon}^p(r)$')
        ax4.legend(loc='best', fontsize=8.5)
        ax4.grid(True, linestyle='--', alpha=0.5)

        fig.suptitle(
            f'2D Elastoplastic BEM Validation: Thick-Walled Cylinder ($P = {pressure:.1f}$ Pa, $\\sigma_Y = {sigma_y:.1f}$ Pa)',
            fontsize=13.5,
        )
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

        all_uy_num = np.concatenate([dpbou[:, 1], dpint[:, 1]])
        all_uy_ex = np.concatenate([exact_bou['uy'], exact_int['uy']])
        err_uy = np.linalg.norm(all_uy_num - all_uy_ex) / np.linalg.norm(all_uy_ex)

        all_vm_num = np.concatenate([stbou[:, 3], stint[:, 3]])
        all_vm_ex = np.concatenate([exact_bou['sigma_vm'], exact_int['sigma_vm']])
        err_vm = np.linalg.norm(all_vm_num - all_vm_ex) / np.linalg.norm(all_vm_ex)

        print("=" * 78)
        print("ELQUABE 2D VALIDATION SUMMARY (PLATE UNDER UNIAXIAL TRACTION)")
        print("=" * 78)
        print(f"Relative L2 Error in Vertical Displacement u_y      : {err_uy:.4e}")
        print(f"Relative L2 Error in von Mises Stress (All Nodes)   : {err_vm:.4e}")
        print("=" * 78)

        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        ax1.plot(y_cont, exact_cont['uy'], 'r-', linewidth=2, label='Analytical Solution $u_y(y)$')
        ax1.plot(bou_coords[:, 1], dpbou[:, 1], 'bo', markersize=7, label='BEM Boundary Nodes')
        ax1.plot(int_coords[:, 1], dpint[:, 1], 'gs', markersize=7, label='BEM Internal Points')
        ax1.set_xlabel('Vertical Coordinate $y$ (m)')
        ax1.set_ylabel('Vertical Displacement $u_y$ (m)')
        ax1.set_title('Vertical Displacement $u_y$ vs. $y$')
        ax1.legend(loc='best')
        ax1.grid(True, linestyle='--', alpha=0.5)

        ax2.plot(y_cont, exact_cont['sigma_vm'], 'r-', linewidth=2, label=f'Analytical $\\sigma_{{\\mathrm{{vM}}}} = \\sigma_y = {traction_y:.1f}$')
        ax2.plot(y_cont, exact_cont['sigma_x'], 'k--', linewidth=1.5, label='Analytical $\\sigma_x = \\tau_{xy} = 0.0$')
        ax2.plot(bou_coords[:, 1], stbou[:, 3], 'bo', markersize=7, label='BEM Boundary $\\sigma_{\\mathrm{vM}}$')
        ax2.plot(int_coords[:, 1], stint[:, 3], 'gs', markersize=7, label='BEM Internal $\\sigma_{\\mathrm{vM}}$')
        ax2.set_xlabel('Vertical Coordinate $y$ (m)')
        ax2.set_ylabel('Stress (Pa)')
        ax2.set_title('Boundary & Internal Stresses vs. $y$')
        ax2.legend(loc='best', fontsize=8.5)
        ax2.grid(True, linestyle='--', alpha=0.5)

        fig.suptitle('ELQUABE 2D Validation: Analytical vs. BEM', fontsize=13)
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
    epstn_png_path='plastic_strain_map.png',
):
    """
    Maps the BEM boundary and internal displacements, stresses, and equivalent
    plastic strain onto the Gmsh quadratic cell mesh (`triangle6` or `quad8`),
    exports a VTK file, and plots the 2D contour maps for displacement,
    von Mises stress, and equivalent plastic strain.
    """
    if inp_params is None:
        inp_params = input_data()
    computed_data = compute_inodes(inp_params['file_name'], inp_params['bound_cond'])

    coord = computed_data['coordinates']
    ordered_boundary_nodes = computed_data['ordered_boundary_nodes']
    cells_conn = computed_data['domain_cells']
    nodc = computed_data['nodc']
    int_nodes = computed_data['internal_nodes']

    n_total = coord.shape[0]
    disp_vec = np.zeros((n_total, 3), dtype=np.float64)
    umag_field = np.zeros(n_total, dtype=np.float64)
    uy_field = np.zeros(n_total, dtype=np.float64)
    sx_field = np.zeros(n_total, dtype=np.float64)
    txy_field = np.zeros(n_total, dtype=np.float64)
    sy_field = np.zeros(n_total, dtype=np.float64)
    vm_field = np.zeros(n_total, dtype=np.float64)
    epstn_field = np.zeros(n_total, dtype=np.float64)

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
        if stbou.shape[1] >= 5:
            epstn_field[node_id] = stbou[bou_idx, 4]

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
        if stint.shape[1] >= 5:
            epstn_field[node_id] = stint[k, 4]

    # 1. Export VTK with quadratic cells (triangle6 or quad8)
    cell_type_name = 'triangle6' if nodc == 6 else 'quad8'
    cells = [(cell_type_name, cells_conn)]
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
            'Equivalent_Plastic_Strain': epstn_field,
        },
    )
    meshio.write(vtk_path, mesh_out)
    print(f"VTK file '{vtk_path}' saved successfully.")

    # 2. Subdivide each quadratic cell into linear triangles for Matplotlib contour plots
    tri3_list = []
    if nodc == 6:
        for t in cells_conn:
            v0, v1, v2, v3, v4, v5 = t
            tri3_list.append([v0, v3, v5])
            tri3_list.append([v3, v1, v4])
            tri3_list.append([v5, v4, v2])
            tri3_list.append([v3, v4, v5])
    else:
        for q in cells_conn:
            v0, v1, v2, v3, v4, v5, v6, v7 = q
            tri3_list.append([v0, v4, v7])
            tri3_list.append([v4, v1, v5])
            tri3_list.append([v5, v2, v6])
            tri3_list.append([v7, v6, v3])
            tri3_list.append([v4, v5, v7])
            tri3_list.append([v5, v6, v7])

    tri3_array = np.array(tri3_list)
    triang = mtri.Triangulation(coord[:, 0], coord[:, 1], tri3_array)

    problem_type = inp_params.get('problem_type', 'cylinder')
    c_radius = None
    if problem_type == 'cylinder':
        params = inp_params.get('analytical_params', {})
        a = params.get('a', 1.0)
        b = params.get('b', 1.5)
        pressure = params.get('pressure', 10.0)
        sigma_y = inp_params.get('plasticity', {}).get('uniax', 24.0)
        ep_sol = compute_cylinder_elastoplastic_solution(
            np.array([[a, 0.0]]), inp_params['E'], inp_params['nu'], sigma_y, a=a, b=b, pressure=pressure
        )
        c_radius = ep_sol['c_radius']
        theta_arc = np.linspace(0.0, 0.5 * np.pi, 100)
        x_arc = c_radius * np.cos(theta_arc)
        y_arc = c_radius * np.sin(theta_arc)

    # 2a. Plot Displacement Contour Map
    if problem_type == 'cylinder':
        disp_plot_field = umag_field
        disp_label = 'Radial Displacement $u_r$ (m)'
        disp_title = 'Radial Displacement $u_r$ (2D Elastoplastic BEM)'
    else:
        disp_plot_field = uy_field
        disp_label = 'Vertical Displacement $u_y$ (m)'
        disp_title = 'Vertical Displacement $u_y$'

    fig1, ax1 = plt.subplots(figsize=(6.5, 5.5))
    cf1 = ax1.tricontourf(triang, disp_plot_field, levels=20, cmap=cm.jet)
    ax1.triplot(triang, color=(0.0, 0.0, 0.0), linewidth=0.35, alpha=0.6)
    if c_radius is not None:
        ax1.plot(x_arc, y_arc, 'w--', linewidth=2.0, label=f'Analytical $c = {c_radius:.3f}$ m')
        ax1.legend(loc='upper right', fontsize=8.5)
    fig1.colorbar(cf1, ax=ax1, label=disp_label)
    ax1.set_title(disp_title)
    ax1.set_xlabel('$x$ (m)')
    ax1.set_ylabel('$y$ (m)')
    ax1.set_aspect('equal')
    plt.savefig(disp_png_path, dpi=300, bbox_inches='tight')
    print(f"Displacement contour map '{disp_png_path}' saved successfully.")
    plt.close(fig1)

    # 2b. Plot von Mises Equivalent Stress Contour Map with Mesh Overlay
    fig2, ax2 = plt.subplots(figsize=(6.5, 5.5))
    cf2 = ax2.tricontourf(triang, vm_field, levels=20, cmap=cm.jet)
    ax2.triplot(triang, color=(0.0, 0.0, 0.0), linewidth=0.35, alpha=0.6)
    ax2.plot(coord[:, 0], coord[:, 1], 'ko', markersize=2.2)
    if c_radius is not None:
        ax2.plot(x_arc, y_arc, 'w--', linewidth=2.0, label=f'Analytical $c = {c_radius:.3f}$ m')
        ax2.legend(loc='upper right', fontsize=8.5)
    fig2.colorbar(cf2, ax=ax2, label='von Mises Stress $\\sigma_{\\mathrm{vM}}$ (Pa)')
    ax2.set_title('von Mises Equivalent Stress $\\sigma_{\\mathrm{vM}}$')
    ax2.set_xlabel('$x$ (m)')
    ax2.set_ylabel('$y$ (m)')
    ax2.set_aspect('equal')
    plt.savefig(vm_png_path, dpi=300, bbox_inches='tight')
    print(f"von Mises contour map '{vm_png_path}' saved successfully.")
    plt.close(fig2)

    # 2c. Plot Equivalent Plastic Strain Contour Map with Mesh Overlay
    fig3, ax3 = plt.subplots(figsize=(6.5, 5.5))
    cf3 = ax3.tricontourf(triang, epstn_field, levels=20, cmap=cm.plasma)
    ax3.triplot(triang, color=(0.2, 0.2, 0.2), linewidth=0.35, alpha=0.5)
    ax3.plot(coord[:, 0], coord[:, 1], 'ko', markersize=2.2)
    if c_radius is not None:
        ax3.plot(x_arc, y_arc, 'c--', linewidth=2.0, label=f'Analytical $c = {c_radius:.3f}$ m')
        ax3.legend(loc='upper right', fontsize=8.5)
    fig3.colorbar(cf3, ax=ax3, label='Equivalent Plastic Strain $\\bar{\\varepsilon}^p$')
    ax3.set_title('Equivalent Plastic Strain $\\bar{\\varepsilon}^p$ Distribution')
    ax3.set_xlabel('$x$ (m)')
    ax3.set_ylabel('$y$ (m)')
    ax3.set_aspect('equal')
    plt.savefig(epstn_png_path, dpi=300, bbox_inches='tight')
    print(f"Plastic strain contour map '{epstn_png_path}' saved successfully.")
    plt.close(fig3)


if __name__ == '__main__':
    inp_params = input_data()
    inp = read_elquabe_input('eq12cb.dat')
    out = read_elquabe_outputs()
    plot_validation(inp, out, inp_params=inp_params, save_path='result_plot.png')
    msh_name = inp_params['file_name']
    if os.path.exists(msh_name + '.msh') or os.path.exists(msh_name + '.geo'):
        export_vtk_and_contours(out, inp_params=inp_params)
