def input_data(problem=1):
    """
    Defines input parameters for the 2D Elastoplastic Boundary Element Method (BEM)
    simulation using continuous quadratic boundary elements and quadratic internal
    cells (6-node triangles or 8-node quadrilaterals).

    Available problems:
        1 : Quarter thick-walled cylinder under internal pressure (placa.geo / placa.msh)
        2 : Square plate under uniaxial traction (placa2.geo / placa2.msh)

    Boundary condition types in `bound_cond`:
        * type_x, type_y = 0 : Prescribed displacement in global x or y direction
        * type_x, type_y = 1 : Prescribed traction in global x or y direction
        * type_x = 2         : Prescribed normal traction/pressure (value_x = t_n);
                               decomposed into global (t_x, t_y) components in preproc.py

    Plasticity parameters in `plasticity`:
        * nsig  : Number of stress components (4 = Plane Strain, 3 = Plane Stress)
        * tolgp : Gauss integration tolerance for internal cells
        * nauto : Gauss integration control (<0: adaptive sub-cells with |nauto| singular points)
        * ncrit : Yield criterion (1 = Tresca, 2 = von Mises, 3 = Mohr-Coulomb, 4 = Drucker-Prager)
        * uniax : Uniaxial yield stress (sigma_Y)
        * frict : Friction angle in degrees (for Mohr-Coulomb / Drucker-Prager)
        * hardi : Isotropic plastic hardening modulus (H')
        * hardk : Kinematic plastic hardening modulus
        * em    : Mixed hardening parameter (1.0 = isotropic, 0.0 = kinematic)
        * nincs : Number of load increments (0 = automatic default from BEM solver)
        * fincr : Load increment scaling factor
        * miter : Maximum iterations per load increment
        * toler : Convergence tolerance for plastic multiplier iterations
    """
    # Select which problem to run (1 = Cylinder under internal pressure, 2 = Plate under traction)
    selected_problem = problem

    if selected_problem == 1:
        # Problem 1: One-quarter thick cylinder under internal pressure (Elastoplastic)
        # Symmetry at y = 0 ('fix1'): t_x = 0 (type_x = 1), u_y = 0 (type_y = 0)
        # Symmetry at x = 0 ('fix2'): u_x = 0 (type_x = 0), t_y = 0 (type_y = 1)
        # Outer radius r = b ('free'): traction-free (t_x = 0, t_y = 0)
        # Inner radius r = a ('loaded'): normal internal pressure (type_x = 2, value_x = -10.0)
        bound_cond = {
            'fix1':   {'type_x': 1, 'value_x': 0.0,   'type_y': 0, 'value_y': 0.0},
            'fix2':   {'type_x': 0, 'value_x': 0.0,   'type_y': 1, 'value_y': 0.0},
            'free':   {'type_x': 1, 'value_x': 0.0,   'type_y': 1, 'value_y': 0.0},
            'loaded': {'type_x': 2, 'value_x': -10.0, 'type_y': 0, 'value_y': 0.0},
        }
        E = 70.0e9
        nu = 0.3
        file_name = 'placa'
        title = 'ELASTOPLASTIC CYLINDER UNDER INTERNAL PRESSURE'
        problem_type = 'cylinder'
        plasticity = {
            'nsig':  4,        # Plane strain
            'tolgp': 1.0e-5,   # Gauss integration tolerance
            'nauto': -6,       # Adaptive cell integration (6 Gauss points for singular integrals)
            'ncrit': 2,        # von Mises yield criterion
            'uniax': 24.0,     # Yield stress sigma_Y (elastic max sigma_eq at r=a is ~31.34)
            'frict': 0.0,      # Friction angle (degrees)
            'hardi': 0.0,      # Perfectly plastic (H' = 0)
            'hardk': 0.0,      # Kinematic hardening modulus
            'em':    1.0,      # Purely isotropic factor
            'nincs': 0,        # 0 = automatic number of increments
            'fincr': 0.95,     # Increment progression factor
            'miter': 50,       # Max plastic iterations per increment
            'toler': 1.0e-5,   # Plastic convergence tolerance
        }
        analytical_params = {
            'a': 1.0,
            'b': 1.5,
            'pressure': 10.0,
            'yield_stress': plasticity['uniax'],
            'ncrit': plasticity['ncrit'],
            'hardi': plasticity['hardi'],
        }

    elif selected_problem == 2:
        # Problem 2: Square plate under uniaxial vertical traction
        bound_cond = {
            'botton': {'type_x': 0, 'value_x': 0.0, 'type_y': 0, 'value_y': 0.0},
            'right':  {'type_x': 1, 'value_x': 0.0, 'type_y': 1, 'value_y': 0.0},
            'top':    {'type_x': 1, 'value_x': 0.0, 'type_y': 1, 'value_y': 1.0},
            'left':   {'type_x': 1, 'value_x': 0.0, 'type_y': 1, 'value_y': 0.0},
        }
        E = 1.0
        nu = 0.0
        file_name = 'placa2'
        title = 'SQUARE PLATE UNDER UNIAXIAL TRACTION'
        problem_type = 'plate'
        plasticity = {
            'nsig':  4,
            'tolgp': 1.0e-5,
            'nauto': -6,
            'ncrit': 2,
            'uniax': 0.8,
            'frict': 0.0,
            'hardi': 0.2,
            'hardk': 0.0,
            'em':    1.0,
            'nincs': 0,
            'fincr': 0.95,
            'miter': 50,
            'toler': 1.0e-5,
        }
        analytical_params = {
            'traction_y': 1.0,
            'yield_stress': plasticity['uniax'],
            'hardi': plasticity['hardi'],
        }

    else:
        raise ValueError(f"Unknown problem selection: {selected_problem}. Choose 1 or 2.")

    return {
        'problem': selected_problem,
        'problem_type': problem_type,
        'title': title,
        'bound_cond': bound_cond,
        'E': E,
        'nu': nu,
        'file_name': file_name,
        'plasticity': plasticity,
        'analytical_params': analytical_params,
    }
