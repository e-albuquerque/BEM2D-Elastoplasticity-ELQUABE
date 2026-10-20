import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection

import meshio

def input_data():
    """
    This function defines input parameters for the Boundary Element Method (BEM)
    simulation applied to the plane elasticity equation with radial integration for domain
    integrals.

    Returns:
        dict: A dictionary containing the following key-value pairs:
            * bound_cond (dict): Boundary conditions for each edge:
                *
            * E (float): Material elastic modulus
            * nu (float): Poisson's ratio
            * file_name (str): Name of the .msh file (Gmsh mesh)
            * qpoint (str): String representing the coordinates of the heat source point
    """
    # Directly define values here for scripting (modify as needed)
    # Boundary conditions are specified in local coordinate system (normal-tangent)
    # type_x = type of the boundary condition in x direction
    # type_y = type of the boundary condition in y direction
    #    (type_x and type_y: 0 = displacement is know while 1 = traction is known)

    bound_cond = {'botton': {'type_x': 0, 'value_x': 0, 'type_y': 0, 'value_y': 0},
                  'right': {'type_x': 1, 'value_x': 0., 'type_y': 1, 'value_y': 0},
        'top': {'type_x': 1, 'value_x': 0, 'type_y': 1, 'value_y': 1},
        'left': {'type_x': 1, 'value_x': 0, 'type_y': 1, 'value_y': 0}
    }
    E = 1. # # Elastic modulus
    nu = 0. # Poisson ration
    file_name = 'placa'


    return {'bound_cond': bound_cond, 'E': E, 'nu': nu, 'file_name': file_name}
def compute_inodes(file_name, bound_cond):
    """
    This function reads a mesh file and extracts information about nodes, elements,
    and boundary conditions.

    Args:
      file_name: Name of the mesh file.
      bound_cond: Dictionary containing boundary conditions for each edge.

    Returns:
      A dictionary containing information about nodes, elements, and boundary
      conditions.
    """
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
            'edge' : edges[key]
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
        'bc_info': bc_info
    }

def mount_bcs(segments, bc_info):
    """
    Constructs a boundary condition matrix from segment information and boundary condition data.

    This function takes a list of segment indices (`segments`) and a dictionary of boundary condition
    information (`bc_info`) to create a matrix (`bcs`) that stores the type and value of the boundary
    condition for each element.

    Args:
        segments (list or np.ndarray): A list or array of integers representing the indices of the boundary
                                       elements (segments) in the problem.
        bc_info (dict): A dictionary containing boundary condition information. The keys are typically
                        strings representing boundary condition names, and the values are dictionaries
                        with the following structure:
                        {
                            'type_x': bc_type,     # The type of boundary condition in x (0 for Dirichlet, 1 for Neumann)
                            'value_x': bc_value    # The value of the boundary condition in x (displacement or traction)
                            'segment': segment_index,  # The index of the segment this BC applies to
                            'type_y': bc_type,          # The type of boundary condition in y (0 for Dirichlet, 1 for Neumann)
                            'value_y': bc_value         # The value of the boundary condition in y (displacement or traction)
                        }

    Returns:
        bcs (np.ndarray): A 2D NumPy array where each row corresponds to a boundary element. The first
                         and third column contains the boundary condition type (0 or 1), and the second
                         and fourth column contains the boundary condition value.
    """

    num_elements = len(segments)

    # Initialize boundary condition matrix
    type_bcs = np.zeros((num_elements, 6), dtype=np.int64)
    val_bcs = np.zeros((num_elements, 6), dtype=np.float64)

    # Populate boundary condition matrix (bcs)
    for i, segment_index in enumerate(segments):
        # Find the boundary condition data corresponding to the current segment
        for bc_name, bc_data in bc_info.items():
            if bc_data['segment'] == segment_index:
                tp_x=bc_data['type_x']
                vl_x= bc_data['value_x']
                tp_y=bc_data['type_y']
                vl_y= bc_data['value_y']
                type_bcs[i,:] = np.array([tp_x,tp_y,tp_x,tp_y,tp_x,tp_y])
                val_bcs[i,:] = np.array([vl_x,vl_y,vl_x,vl_y,vl_x,vl_y])
                break  # Move on to the next segment once the BC is found
    return type_bcs,val_bcs

def comp_xyplot(elem, nodes):
    """
    Calculates the nodes and outward normal vector for each boundary element.

    This function takes element connectivity information (`elem`) and node coordinates (`nodes`)
    to determine the geometric center of each boundary element and the unit normal vector
    pointing outward from that element.

    Args:
        elem (np.ndarray): A 2D array defining the connectivity of the boundary elements.
                          Each row contains the indices of the two nodes that form an element.
                          Shape: (num_elements, 2).
        nodes (np.ndarray): A 2D array containing the coordinates (x, y) of each node in the mesh.
                           Shape: (num_nodes, 2).

    Returns:
        node_des (np.ndarray): A 2D array containing the coordinates (x, y) of the discontinuous nodes
                              of each boundary element. Shape: (2*num_elements, 2).
        normal (np.ndarray): A 2D array containing the unit outward normal vectors for each
                             boundary element. Shape: (num_elements, 2).
    """
    npp = 8                              # Number of plot points
    eta = np.linspace(-1., 1., npp)        # Plot points

    num_elements = elem.shape[0]

    xplot = np.zeros((num_elements,npp))        # Store x-coordinates for plotting
    yplot = np.zeros((num_elements,npp))        # Store y-coordinates for plotting


      # Define the local coordinates (xi) for nodes within a quadratic element
    xi = np.array([-1., 0., 1.])  # Gauss-Lobatto points for quadratic elements

    # Iterate over each boundary element
    for i in range(num_elements):
        # Get node indices for the current element
        node1, node3, node2 = elem[i]

        # Get coordinates of the two nodes forming the element
        x1, y1 = nodes[node1][0:2]
        x2, y2 = nodes[node2][0:2]
        x3, y3 = nodes[node3][0:2]

        # Iterate over Gauss points for numerical integration
        for k in range(npp):
            # Shape functions and their derivatives
            N1 = 0.5 * eta[k] * (eta[k] - 1)
            N2 = 1 - eta[k]**2
            N3 = 0.5 * eta[k] * (eta[k] + 1)

            # Calculate coordinates and normal vector at the Gauss point
            x, y = N1 * x1 + N2 * x2 + N3 * x3, N1 * y1 + N2 * y2 + N3 * y3
            ivet = npp * i + k  # Index for storing plot coordinates
            xplot[i,k], yplot[i,k] = x, y

    return xplot,yplot
def show_geometry(nodes,xplot,yplot,elem):
    xnode=nodes[:,0]
    ynode=nodes[:,1]
    ax = plt.axes()
    plt.plot(xnode[elem[:,0]],ynode[elem[:,0]],"gd",markersize=8)	# Plot the node of the elements
    plt.plot(xnode,ynode,"ro",markersize=2)	# Plot the node of the elements
    plt.axis("equal")

    n_el = elem.shape[0]
    for i in range(n_el):
        x=xplot[i,:]
        y=yplot[i,:]
        line_segments = LineCollection([list(zip(x, y))],
                                        linestyles='solid')
        ax.add_collection(line_segments)

    plt.show()

def format_coord(val):
    # Arredonda para evitar problemas de preciso numrica do Gmsh (ex: 2.999999999999 -> 3.0)
    val_rounded = round(val, 4)
    if val_rounded.is_integer():
        return f"{int(val_rounded)}."
    return f"{val_rounded}"

def create_input_file(filename, computed_data, int_points, \
                      type_bcs, val_bcs, input_params):
    """
    Creates an input file for the Fortran BEM program based on computed mesh data
    and boundary conditions.
    """
    with open(filename, 'w') as f:
        # 1. Title Line
        f.write("SQUARE PLATE UNDER TIME-HARMONIC TRACTION (4 QUADRATIC ELEMENTS)\n")

        # Extract data from computed_data and input_params
        nodes = computed_data['coordinates']
        elements = computed_data['elements']
        num_elements = elements.shape[0]
        num_internal_points = len(int_points)


        shear_modulus_real = input_params['E'] / (2 * (1 + input_params['nu']))

        poisson_ratio = input_params['nu']

        # 2. Basic Parameter Line (formatted to match squa4.dat)
        f.write(f" {num_elements},{num_internal_points},\
        {shear_modulus_real:.1f},{poisson_ratio:.2f}\n")

        # 3. Boundary Nodes Coordinates Lines
        # Sequência geométrica desejada: no1 e no_meio de cada elemento
        node_coords_list = []
        for elem in elements:
            no1_idx = elem[0]
            no_meio_idx = elem[2]
            x1, y1 = nodes[no1_idx][0], nodes[no1_idx][1]
            xm, ym = nodes[no_meio_idx][0], nodes[no_meio_idx][1]
            node_coords_list.append(f"{format_coord(x1)}  {format_coord(y1)}")
            node_coords_list.append(f"{format_coord(xm)}  {format_coord(ym)}")

        # Junta os pares usando o mesmo espaçamento largo do arquivo squa4.dat original
        node_coords_str = "  ".join(node_coords_list)
        f.write("  " + node_coords_str + "\n")

        # 5. Boundary Conditions Lines
        for i in range(num_elements):
            # X and Y direction BCs for 3 nodes
            f.write(f"{type_bcs[i, 0]} {val_bcs[i, 0]:.1f} ")
            f.write(f"{type_bcs[i, 1]} {val_bcs[i, 1]:.1f} ")
            f.write(f"{type_bcs[i, 2]} {val_bcs[i, 2]:.1f} ")
            f.write(f"{type_bcs[i, 3]} {val_bcs[i, 3]:.1f} ")
            f.write(f"{type_bcs[i, 4]} {val_bcs[i, 4]:.1f} ")
            f.write(f"{type_bcs[i, 5]} {val_bcs[i, 5]:.1f} \n")

        # 6. Internal Points Coordinates Lines
        int_points_str = " ".join(f"{nodes[p][0]:.1f}  {nodes[p][1]:.1f}" for p in int_points)
        f.write(" " + int_points_str + "\n")



if __name__ == '__main__':
    # 1. Executa pré-processamento para gerar o arquivo cube8.dat (nó 770 rastreado)
    # Read information about input data
    inp_data = input_data()
    # Format input data
    computed_data = compute_inodes(inp_data['file_name'], inp_data['bound_cond'])
    # Generate bcs array with boundary conditions in each element
    type_bcs,val_bcs = mount_bcs(computed_data['segments'], computed_data['bc_info'])
    xplot,yplot = comp_xyplot(computed_data['elements'], computed_data['coordinates'])
    show_geometry(computed_data['coordinates'],xplot,yplot,computed_data['elements'],)
    # Re-get the global `inp_data` dictionary to avoid shadowing
    global_input_data = input_data()
    inp_data_filename='eq12cb.dat'
    int_points=computed_data['internal_nodes']
    # Call the function with the correct arguments
    create_input_file(inp_data_filename, computed_data, int_points, \
                      type_bcs, val_bcs, global_input_data)