"""
Created on Mon Feb 19 09:46:43 2023

@author: lzampa
"""
# -----------------------------------------------------------------------------
## Importing base python modules

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

# Close all open figures
plt.close('all') 

# -----------------------------------------------------------------------------
def line_plane_intersection(plane_point, plane_normal, line_point1, line_point2):
    """
    Calculate the intersection point between a line and a plane in 3D space.

    Parameters:
    plane_point (np.array): A point on the plane, represented as a numpy array of three elements.
    plane_normal (np.array): The normal vector to the plane, represented as a numpy array of three elements.
    line_point1, line_point2 (np.array): Two points on the line, each represented as a numpy array of three elements.

    Returns:
    intersection_point (np.array): The intersection point between the line and the plane, if it exists, otherwise None.
    """

    # Convert inputs to numpy arrays for easier vector operations
    plane_normal = np.asarray(plane_normal)
    plane_point = np.asarray(plane_point)
    line_point1 = np.asarray(line_point1)
    line_point2 = np.asarray(line_point2)
    
    # Calculate the vector of the line
    line_vector = line_point2 - line_point1
    # Calculate the vector from the point on the line to the point on the plane
    plane_vector = plane_point - line_point1
    
    # Calculate the parameter t that represents the intersection point along the line
    t = np.dot(plane_normal, plane_vector) / np.dot(plane_normal, line_vector)

    # If t is between 0 and 1, the intersection point is on the line between line_point1 and line_point2
    if 0 <= t <= 1:
        # Calculate the intersection point
        intersection_point = line_point1 + t * line_vector
        return intersection_point
    else:
        # If t is not between 0 and 1, there is no intersection between the line and the plane
        return None

# -----------------------------------------------------------------------------
def project_point_onto_plane(point, plane_point, plane_normal):
    """
    Project a point onto a plane in 3D space.

    Parameters:
    point (np.array): The point to be projected, represented as a numpy array of three elements.
    plane_point (np.array): A point on the plane, represented as a numpy array of three elements.
    plane_normal (np.array): The normal vector to the plane, represented as a numpy array of three elements.

    Returns:
    projected_point (np.array): The projection of the input point onto the plane.
    """

    # Convert inputs to numpy arrays for easier vector operations
    point = np.asarray(point)
    plane_point = np.asarray(plane_point)
    plane_normal = np.asarray(plane_normal)

    # Calculate the vector from the point on the plane to the input point
    point_vector = point - plane_point
    # Calculate the projection of point_vector onto the plane
    projection_vector = point_vector - np.dot(point_vector, plane_normal) * plane_normal / np.linalg.norm(plane_normal)**2

    # Calculate the projected point
    projected_point = plane_point + projection_vector

    return projected_point

# -----------------------------------------------------------------------------
def min_travel_time_SinglePlane( V1, V2, source_xyz, receiver_xyz, plane_point, 
                                 plane_normal, space_threshold=1e-6 ):
    """
    Calculate the minimum travel time for a seismic wave to travel from a source to a receiver, 
    refracting at an interface defined by a plane.

    Parameters:
    V1, V2 (float): The seismic velocities above and below the interface, respectively.
    source_xyz, receiver_xyz (np.array): The coordinates of the source and receiver, respectively.
    plane_point (np.array): A point on the interface plane.
    plane_normal (np.array): The normal vector to the interface plane.
    space_threshold (float): The threshold for the spatial accuracy of the solution.

    Returns:
    ray_path (np.array): The coordinates of the ray path from the source, through the interface, to the receiver.
    travel_time (float): The travel time of the seismic wave along the ray path.
    """

    # Convert inputs to numpy arrays for easier vector operations
    source_xyz = np.asarray(source_xyz)
    receiver_xyz = np.asarray(receiver_xyz)
    plane_point = np.asarray(plane_point)

    # We need to check wheather the source or reciever points stand on the plane, 
    # if so, we move the points a little bit to avoid the singularity.
    # To do that we uses the plane equation Ax + By + Cz + D = 0, 
    # where (A, B, C) is the normal vector of the plane, 
    # and D is the distance from the plane to the origin.

    # First, we calculate D in the plane equation
    D = np.dot( plane_normal, plane_point )
    # Then, we substitute the point into the plane equation
    dotp = np.dot(plane_normal, source_xyz) - D
    # Finally, we check if the dot product 
    # of the two parrallel vectors minus D is close to zero
    if np.isclose(dotp, 0) :
        # Point shifted down by 1e-6 in z direction
        source_xyz[2] = source_xyz[2] - 1e-6 

    # Project the source and receiver points onto the plane
    Sp = project_point_onto_plane(source_xyz, plane_point, plane_normal)
    Rp = project_point_onto_plane(receiver_xyz, plane_point, plane_normal)

    # Calculate the intersection point of the line between the source and receiver with the plane
    I = line_plane_intersection(plane_point, plane_normal, source_xyz, receiver_xyz)
    if I is None:
        return None, None

    # Calculate distances for Snell's law and the bisection method
    a = np.linalg.norm(receiver_xyz - Rp)
    b = np.linalg.norm(source_xyz - Sp)
    c = np.linalg.norm(Sp - Rp)

    # Initial values for the bisection method
    x1, x2 = 0.0, c
    x = (x1 + x2) / 2.0

    # Bisection method to find the point of refraction on the interface
    while True:

        i0 = np.copy( I ) 
        
        f1 = x1/(V1 * np.sqrt(a**2 + x1**2)) - (c - x1)/(V2 * np.sqrt(b**2 + (c - x1)**2))
        f2 = x2/(V1 * np.sqrt(a**2 + x2**2)) - (c - x2)/(V2 * np.sqrt(b**2 + (c - x2)**2))
        fx = x/(V1 * np.sqrt(a**2 + x**2)) - (c - x)/(V2 * np.sqrt(b**2 + (c - x)**2) )

        if f1 * fx < 0:
            x2 = x
            x = (x1 + x2) / 2.0
        if f2 * fx < 0 :
            x1 = x
            x = (x1 + x2) / 2.0

        I = Rp + ( Sp - Rp ) * (x/c)

        if np.linalg.norm(i0-I) <= space_threshold :
            break

    # Define the ray path
    x_coords = np.array( [ source_xyz[0], I[0], receiver_xyz[0] ] )
    y_coords = np.array( [ source_xyz[1], I[1], receiver_xyz[1] ] )
    z_coords = np.array( [ source_xyz[2], I[2], receiver_xyz[2] ] )

    ray_path = np.column_stack( ( x_coords, y_coords, z_coords ) )
    travel_time = np.linalg.norm( source_xyz - I ) / V2 + np.linalg.norm( I - receiver_xyz ) / V1

    return ray_path, travel_time

# -----------------------------------------------------------------------------
def min_travel_time_MultiPlane( V, source_xyz, receiver_xyz, plane_points, 
                                plane_normals, time_threshold=1e-6,
                                space_threshold=1e-6, plot=False, prj2plane='xz' ):
    """
    This function calculates the minimum travel time of a seismic wave through multiple planes.

    Parameters:
    V (numpy array): Velocities of the seismic wave in each plane.
    source_xyz (numpy array): Coordinates of the source of the seismic wave.
    receiver_xyz (numpy array): Coordinates of the receiver of the seismic wave.
    plane_points (numpy array): Points on each plane.
    plane_normals (numpy array): Normal vectors to each plane.
    time_threshold (float): Threshold for the difference in travel time for the loop to stop.
    space_threshold (float): Threshold for the difference in space for the loop to stop.
    plot (bool): If True, the function will plot the ray paths.
    prj2plane (str): Projection plane for the plot. Can be 'xz', 'xy', or 'yz'.

    Returns:
    ray_path (numpy array): Path of the seismic wave.
    travel_time (float): Minimum travel time of the seismic wave.
    """

    # Convert inputs to numpy arrays
    V = np.asarray( V )
    plane_points = np.asarray( plane_points )
    plane_normals = np.asarray( plane_normals )
    source_xyz = np.asarray( source_xyz )
    receiver_xyz = np.asarray( receiver_xyz )

    # Make copies of the inputs
    V0 = np.copy( V )
    plane_points0 = np.copy( plane_points )
    plane_normals0 = np.copy( plane_normals )

    # Initialize travel time list and ray path
    travel_time_list = [1e10]
    ray_path = np.full( (V.size+1, 3 ), np.nan )
    ray_path[0,:] = source_xyz
    ray_path[1:,:] = receiver_xyz

    n = 0
    while True :
        
        # Flip ray path after first iteration
        if n!= 0 :
            ray_path = np.flipud( ray_path )

        # Flip velocities, plane points, and plane normals on odd iterations
        if n % 2 == 1 :
            V = np.flipud( V0 )
            plane_points = np.flipud( plane_points0 )
            plane_normals = np.flipud( plane_normals0 )
        else :
            V = np.copy( V0 )
            plane_points = np.copy( plane_points0 )
            plane_normals = np.copy( plane_normals0 )

        travel_time = 0
        for i in range(len(V) - 1):  

            # Get velocity, plane point, and plane normal for current plane
            V2 = V[i]
            V1 = V[i + 1]
            plane_point = plane_points[i]
            plane_normal = plane_normals[i]
            
            # Get source and receiver for current plane
            S = ray_path[i]
            R = ray_path[i+2]

            # Calculate minimum travel time for current plane
            r, t = min_travel_time_SinglePlane(V1, V2, S, R, plane_point, plane_normal, space_threshold)
            if r is None :
                ray_path[i+1,:] = ray_path[i,:]
                continue
            else :
                ray_path[i+1,:] = r[1,: ]
                travel_time += t

        n += 1

        # Add current travel time to list
        travel_time_list.append( travel_time )

        # Break loop if difference in travel time is below threshold and iteration is odd
        if ( np.abs( travel_time_list[-1] - travel_time_list[-2] ) <= time_threshold ) and ( n % 2 == 1) :
            break

    # Remove NaN values from ray path
    inan = np.isnan( ray_path[:,0] )
    ray_path = ray_path[ ~inan ]

    # Plotting 
    # -------
    if plot is True :
        plot_velmod_raypath( source_points, receiver_points, 
                             [ ray_path ], plane_points, prj2plane=prj2plane)
        
    return ray_path, travel_time

# -----------------------------------------------------------------------------
def min_travel_times_paths_1D( V, depths, source_points, receiver_points, 
                               SR_combinations=None, plot=False, prj2plane='xz'):
    """
    This function calculates the minimum travel times and paths for seismic waves in a 1D velocity model.

    Parameters:
    V (array-like): Velocity values.
    depths (array-like): Corresponding depths for the velocity values.
    source_points (array-like): Source points in 3D space (x=1stCol, y=2ndCol, z=3rdCol, name=4thColo[optional]).
    receiver_points (array-like): Receiver points in 3D space (x=1stCol, y=2ndCol, z=3rdCol, name=4thColo[optional]).
    SR_combinations (array-like, optional): Combinations of source and receiver points (Sname=1stCol, Rname=2ndCol). 
    If None, all possible combinations are used.
    plot (bool, optional): If True, the velocity model and ray paths are plotted.
    prj2plane (str, optional): Projection plane for plotting ('xz', 'yz', or 'xy').

    Returns:
    all_ray_paths (list): List of all ray paths.
    all_travel_times (list): List of corresponding travel times.
    """

    # Convert inputs to numpy arrays
    V = np.asarray(V).ravel()
    depths = np.asarray(depths).ravel()
    source_points = np.asarray(source_points)
    receiver_points = np.asarray(receiver_points)

    # Ensure depths and V have the same length
    while len(depths) > len(V):
        V = np.append(V, V[-1])
    while len(depths) < len(V):
        V = V[:-1]

    # Sort depths and V in ascending order of depth
    dep_vel = np.column_stack((depths, V))
    depths = depths[dep_vel[:, 0].argsort()]
    V = V[dep_vel[:, 0].argsort()]

    # Define planes at each depth
    plane_points = []
    plane_normals = []
    for i in range(len(depths) - 1):
        plane_points.append( np.array( ( 0, 0, depths[i] ) ) )
        plane_normals.append( np.array( ( 0, 0, 1 ) ) )

    # Initialize lists to store results
    all_travel_times = []
    all_ray_paths = []

    # Extract or generate source and receiver names
    if source_points.shape[1] == 3:
        source_names = np.arange(0, source_points.shape[0])
    else:
        source_names = source_points[:, 0]
        source_points = source_points[:, :-1]
    if receiver_points.shape[1] == 3:
        receiver_names = np.arange(0, receiver_points.shape[0])
    else:
        receiver_names = receiver_points[:, 0]
        receiver_points = receiver_points[:, :-1]

    # Generate source-receiver combinations if not provided
    if SR_combinations is None:
        SR_combinations = np.dstack(np.meshgrid(source_names, receiver_names)).reshape(-1, 2)

    # Calculate minimum travel time and ray path for each source-receiver pair
    for source_name, receiver_name in SR_combinations:
        source_xyz = source_points[source_names == source_name]
        receiver_xyz = receiver_points[receiver_names == receiver_name]
        ray_path, travel_time = min_travel_time_MultiPlane(V, source_xyz, receiver_xyz, plane_points, plane_normals)
        ray_path  = np.column_stack( ( ray_path, 
                                       np.full( (ray_path.shape[0]), source_name ),  
                                       np.full( (ray_path.shape[0]), receiver_name ) ) )
        travel_time  = np.column_stack( ( travel_time, source_name, receiver_name ) ) 
        all_ray_paths.append(ray_path)
        all_travel_times.append(travel_time)

    # Plot results if requested
    if plot is True:
        plot_velmod_raypath(source_points, receiver_points, plane_points, all_ray_paths, prj2plane=prj2plane)

    return all_ray_paths, all_travel_times

# --------------------------------------------------------------------------------
def plot_velmod_raypath( source_points, receiver_points, 
                         plane_points, ray_paths=None, prj2plane='xz') :
    """
    This function plots the velocity model and ray paths of a seismic wave in a 1D medium.

    Parameters:
    source_points (list): List of source points of the seismic wave.
    receiver_points (list): List of receiver points of the seismic wave.
    plane_points (list): List of points on each plane.
    prj2plane (str): Projection plane for the plot. Can be 'xz', 'xy', or 'yz'.

    Returns:
    None
    """

    # Close all previous plots
    plt.close('all')

    # Create new figure and axes
    fig, ax = plt.subplots()

    # Get colormap and normalize it to the range of velocities
    cmap = plt.get_cmap('viridis')
    norm = mcolors.Normalize(vmin=min(V), vmax=max(V))

    # Set projection plane
    if prj2plane == 'xz':
        n1 = 0
        n2 = 2
    elif prj2plane == 'xy':
        n1 = 0
        n2 = 1
    elif prj2plane == 'yz':
        n1 = 1
        n2 = 2

    # Get all points and their minimum and maximum z-coordinates
    all_points = np.vstack( ( source_points, plane_points, receiver_points ) )
    min_z, max_z = np.min( all_points[:,n2] ), np.max( all_points[:,n2] )
    addz = ( max_z - min_z ) * 0.1

    # Sort plane points and add padding
    plane_points2 = np.sort(np.array( plane_points )[:,n2], axis=0)
    plane_points2 = np.insert(plane_points2, 0, min_z-addz) 
    plane_points2 = np.append(plane_points2, max_z+addz) 

    # Plot ray paths
    unique_labels = set()
    for ray_path in ray_paths:
        plt.plot( ray_path[:,n1], ray_path[:,n2], linestyle='-', color='b', 
                  label='ray paths' if 'ray paths' not in unique_labels else "", zorder=2)
        unique_labels.add('ray paths')

    # Plot velocity model
    xarr = np.linspace(ax.get_xlim()[0], ax.get_xlim()[1], 100)
    for i in range( len(V) ) :
        color = cmap(norm(V[i]))
        plt.fill_between(xarr, [plane_points2[i]], [plane_points2[i+1]], color=color, alpha=0.5 )

    # Add colorbar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cmap_transparent = mcolors.LinearSegmentedColormap.from_list(
        'transparent_viridis', 
        [(r, g, b, 0.5) for r, g, b, _ in cmap(np.linspace(0, 1, cmap.N))])
    sm_transparent = plt.cm.ScalarMappable(cmap=cmap_transparent, norm=norm)
    sm_transparent.set_array([])
    fig.colorbar(sm_transparent, ax=ax, label='Velocity')

    # Plot source and receiver points
    for i, source_point in enumerate( source_points ):
        plt.scatter( source_point[n1], source_point[n2], c='r', marker='*', 
                     label='sources' if 'sources' not in unique_labels else "", zorder=3 )
        unique_labels.add('sources')
        if np.size( source_point ) > 3 :
            plt.text( source_point[n1], source_point[n2], str( source_point[3] ), zorder=5 )
        else :
            plt.text( source_point[n1], source_point[n2], str( i ), zorder=5 )
    
    for i, receiver_point in enumerate( receiver_points ):
        plt.scatter( receiver_point[n1], receiver_point[n2], c='r', marker='v', 
                     label='recievers' if 'recievers' not in unique_labels else "", zorder=3 )
        unique_labels.add('recievers')
        if np.size( receiver_point ) > 3 :
            plt.text( receiver_point[n1], receiver_point[n2], str( receiver_point[3] ), zorder=5 )
        else :
            plt.text( receiver_point[n1], receiver_point[n2], str( i ), zorder=5 )

    # Add legend, grid, and adjust layout
    plt.legend()
    plt.grid()
    plt.tight_layout()

    # Show plot
    plt.show()
    
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------
# --------------------------------------------------------------------------------

# Main program
if __name__ == "__main__":

    # Define velocities (V), source/receiver positions (x, z), and interface depths
    V = [ 8000, 5000, 3000, 2000 ]
    depths = [ -4000, -2000, -1000, 0 ]
    source_points = [ ( -1000, 0, -5000 ), ( -2000, 0, -3500 ) ]  
    receiver_points = [ ( 1000, 0, 0 ), ( 500, 0, 0 ) ]  

    # Calculate the minimum travel time path
    all_ray_paths, all_travel_times = min_travel_times_paths_1D( V, depths,
                                source_points, receiver_points, plot=True, prj2plane='xz' )
    

