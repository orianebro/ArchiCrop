import numpy as np
from openalea.mtg import MTG

## Initialize parameters

# Fixed parameters 
reproductive_start_time = 50  # Example time for reproductive organ appearance
max_leaf_size = 10  # Maximum leaf size constraint

# Parameters we want to vary
phyllochron = 7  # Example phyllochron in days


## Data from crop model ()

# Define functions for plant height and LAI temporal evolution (beta functions)
def height_function(t):
    # Example beta function for height evolution
    return 1.5 * np.exp(-0.02 * t) + 0.5 * np.sin(0.1 * t)

def LAI_function(t):
    # Example beta function for LAI evolution
    return 2 * np.exp(-0.02 * t) + 0.8 * np.cos(0.1 * t)


# Create an empty MTG
plant_mtg = MTG()

# Add main stem as a component to the MTG
main_stem_id = plant_mtg.add_component(plant_mtg.root, label='MainStem')

# Loop over time steps
for t in range(100):  # Example growth period of 100 days
    # Update plant height and LAI
    height = height_function(t)
    LAI = LAI_function(t)
    
    # Check if it's time for a new leaf to appear
    if t % phyllochron == 0:
        # Create a new internode and leaf
        internode_id = plant_mtg.add_child(main_stem_id, label=f'Internode_{t // phyllochron}')
        leaf_id = plant_mtg.add_child(internode_id, label=f'Leaf_{t // phyllochron}')
    
    # Calculate growth demands for each leaf
    for leaf in plant_mtg.components(plant_mtg.root):
        if plant_mtg.label(leaf).startswith('Leaf'):
            leaf_age = t / phyllochron
            leaf_demand = max_leaf_size * (1 - np.exp(-leaf_age))  # Example beta function for growth demand
            # Adjust leaf demand based on age-dependent surface growth demand (beta function)
            # Adjust leaf demand based on LAI increase
            # Distribute available resources among leaves based on their demands
    
    # Update leaf dimensions and angles
    # for leaf in plant_mtg.components(plant_mtg.root):
    #     if plant_mtg.label(leaf).startswith('Leaf'):
            # Calculate angle with axis based on leaf age
            # Update leaf dimensions (length, width) and angles
    
    # Check if it's time for reproductive organs to appear
    # if t >= reproductive_start_time:
        # Create flowers/fruits on appropriate nodes based on certain conditions

# Output final plant structure and characteristics
print("Final plant structure:")
for vid in plant_mtg.vertices():
    print("Vertex:", vid, "Label:", plant_mtg.label(vid))


from openalea.plantgl.all import *

# Create PlantGL shapes for stem (example)
stem_geometry = Cylinder(1.0, 5.0)  # Cylinder(radius, height)
# stem_material = Material(Color3, 0.5, 0.3, 0.1)  # Brown color for stem
PlantGL_shape_for_stem = stem_geometry # + stem_material

# Define PlantGL_shape_for_leaf (example)
# Define vertices and faces for the leaf surface
leaf_vertices = [(0, 0, 0), (1, 0, 0), (0.5, 1, 0)]  # Example vertices
leaf_faces = [(0, 1, 2)]  # Example faces

# Create TriangleSet geometry for the leaf
leaf_geometry = TriangleSet(leaf_vertices, leaf_faces)
# leaf_material = Material(Color3, 0, 0.8, 0)  # Green color for leaf
PlantGL_shape_for_leaf = leaf_geometry # + leaf_material


# Create Shape objects for stem and leaf
stem_shape = Shape(PlantGL_shape_for_stem) 
leaf_shape = Shape(PlantGL_shape_for_leaf) 

# Create materials for stem and leaf (optional)
# stem_material = Material(Color3(0.2, 0.8, 0.2))  # Green color for stem
# leaf_material = Material(Color3(0.8, 0.2, 0.2))  # Red color for leaf
# stem_shape.appearance = stem_material
# leaf_shape.appearance = leaf_material

# Position elements based on their coordinates
stem_shape.setPosition(0, 0, 0)  # Example position for stem
leaf_shape.setPosition(1, 0, 0)  # Example position for leaf

# Create a scene and add shapes to it
shapes = [stem_shape, leaf_shape]
scene = Scene(shapes)

# Visualize the MTG
Viewer.display((scene))

# # Add components to the MTG (example)
# main_stem_id = plant_mtg.add_component(plant_mtg.root, label='MainStem')
# internode_id = plant_mtg.add_child(main_stem_id, label='Internode')
# leaf_id = plant_mtg.add_child(internode_id, label='Leaf')



# Assign shapes to components (example)
plant_mtg.properties()[main_stem_id] = {'shape': stem_shape}
plant_mtg.properties()[leaf_id] = {'shape': leaf_shape}
