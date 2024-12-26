import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import shapely
import matplotlib.pyplot as plt
import math

import math

import math


def add_and_move_midpoints(nodes, distance):
    """
    Add nodes at the midpoint (based on Euclidean distance) between each pair of consecutive nodes in the input list,
    and then move the new nodes away from the original line by the specified distance along the normal vector.

    Args:
        nodes (list): A list of [x, y] coordinate pairs.
        distance (float): The distance to move the new nodes away from the original line.

    Returns:
        list: A new list of [x, y] coordinate pairs with midpoints added and moved.
    """
    nodes.append(nodes[0])
    new_nodes = []
    for i in range(len(nodes) - 1):
        x1, y1 = nodes[i]
        x2, y2 = nodes[i + 1]

        # Calculate the Euclidean distance between the two nodes
        distance_between = math.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)

        # Calculate the midpoint coordinates
        midpoint_x = x1 + (x2 - x1) / 2
        midpoint_y = y1 + (y2 - y1) / 2

        # Calculate the normal vector
        normal_x = y2 - y1
        normal_y = -(x2 - x1)
        normal_length = math.sqrt(normal_x ** 2 + normal_y ** 2)
        normal_x /= normal_length
        normal_y /= normal_length

        # Move the midpoint along the normal vector
        new_midpoint_x = midpoint_x + normal_x * distance
        new_midpoint_y = midpoint_y + normal_y * distance

        new_nodes.append([x1, y1])
        new_nodes.append([new_midpoint_x, new_midpoint_y])

        # Add the last node
        if i == len(nodes) - 2:
            new_nodes.append([x2, y2])

    return new_nodes


# Create the graph
G = nx.Graph()

# Add nodes
G.add_nodes_from([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])

# Add edges
G.add_edges_from([(1, 2), (2, 3), (3, 4), (4, 5),
                  (6, 7), (7, 8), (8, 9), (9, 10)])

# Define the clusters
cluster1 = [1, 2, 3, 4, 5]
cluster2 = [6, 7, 8, 9, 10]
cluster3 = [1, 2, 3]
cluster4 = [4, 5]
cluster5 = [7, 8, 9]

clusters = [cluster1, cluster2, cluster3, cluster4, cluster5]

# Create a dictionary to store the cluster memberships
cluster_dict = {node: [] for node in G.nodes()}
for node in G.nodes():
    for i, current_cluster in enumerate(clusters):
        cluster_dict[node].append(i+1)

# Define the colors for the clusters
cluster_colors = ['red', 'blue', 'green', 'orange', 'purple']

# Draw the graph
pos = nx.spring_layout(G, k=0.5)


plt.figure(figsize=(10, 10))
# Draw the nodes
node_colors = [cluster_colors[max(cluster_dict[node]) - 1] for node in G.nodes()]
nx.draw_networkx_nodes(G, pos, node_color='gray', edgecolors='black', linewidths=2, alpha=0.7)

# Draw the edges
nx.draw_networkx_edges(G, pos)
# Draw the cluster blobs
for i, cluster in enumerate(clusters):
    cluster_pos = [pos[node] for node in cluster]
    x, y = zip(*cluster_pos)
    hull = list(np.array(list(zip(x, y))))
    print(len(hull))
    hull = add_and_move_midpoints(hull, distance=-0.2)
    print(len(hull))
    print("")
    # Calculate the centroid of the hull
    centroid = np.mean(hull, axis=0)


    # Calculate the scaling factor to extend the hull by 10%
    scaling_factor = 1.1

    # Scale the hull around the centroid
    scaled_hull = [(point - centroid) * scaling_factor + centroid for point in hull]

    polygon = patches.Polygon(scaled_hull, closed=True, fill=True, color=cluster_colors[i], linewidth=0, alpha=0.3)
    plt.gca().add_patch(polygon)

# Add labels
nx.draw_networkx_labels(G, pos)

plt.axis('off')
plt.show()
