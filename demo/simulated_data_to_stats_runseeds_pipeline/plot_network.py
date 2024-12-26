from graphviz import Digraph

dot = Digraph(comment='Group Graph')
dot.graph_attr['compound'] = 'true'

# Define subgraphs (clusters)
with dot.subgraph(name='cluster_0') as c:
    c.attr(color='red')
    c.node_attr.update(style='filled', color='lightgrey')
    c.edges([('1','2'), ('2','3'), ('3','4')])
    c.attr(label='Group A')

with dot.subgraph(name='cluster_1') as c:
    c.attr(color='blue')
    c.node_attr.update(style='filled', color='blue')
    c.edges([('5','6'), ('6','7')])
    c.attr(label='Group B')

# New overlapping cluster_3 with nodes 4 and 7
with dot.subgraph(name='cluster_2') as c:
    c.attr(color='green')
    c.node_attr.update(style='filled', color='green')
    c.node('4')
    c.node('7')
    c.edges([('4', '7')])
    c.attr(label='Group C')

# Separate node 8
dot.node('8')

# Edges from cluster_0 and cluster_1 to node 8
dot.edge('3', '8', ltail='cluster_0')
dot.edge('6', '8', ltail='cluster_1')
dot.edge('4', '9', ltail='cluster_2')

dot.render('group_graph', view=True)