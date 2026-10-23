def input_data(problem=2):
    """
    Defines input parameters for the 2D Elastostatic Boundary Element Method (BEM)
    simulation using continuous quadratic boundary elements (ELQUABE).

    Available problems:
        1 : Quarter cylinder under internal pressure (placa.geo / placa.msh)
        2 : Square plate under uniaxial traction (placa2.geo / placa2.msh)

    Boundary condition types in `bound_cond`:
        * type_x, type_y = 0 : Prescribed displacement in global x or y direction
        * type_x, type_y = 1 : Prescribed traction in global x or y direction
        * type_x = 2         : Prescribed normal traction/pressure (value_x = t_n);
                               decomposed into global (t_x, t_y) components in preproc.py
    """
    # Select which problem to run (1 = Cylinder under internal pressure, 2 = Plate under traction)
    selected_problem = problem

    if selected_problem == 1:
        # Problem 1: One-quarter thick cylinder under internal pressure
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
        title = 'CYLINDER UNDER INTERNAL PRESSURE'
        problem_type = 'cylinder'
        analytical_params = {'a': 1.0, 'b': 1.5, 'pressure': 10.0}

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
        title = 'SQUARE PLATE UNDER UNIAXIAL STATIC TRACTION'
        problem_type = 'plate'
        analytical_params = {'traction_y': 1.0}

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
        'analytical_params': analytical_params,
    }
