# testing of operator-functionality
# for testing to be fully rigurous, the "data"-folder should be deleted before running the tests


import unittest
from WOHairyOperators import WOHairyGC
from sage.all import *
import Shared
import matplotlib.pyplot as plt


def get_partition(n_vertices, n, n_omega, n_epsilon):

    if n_vertices > 0: inner_vertices = [list(range(0, n_vertices))]
    else: inner_vertices = []

    numbered_vertices = [[j] for j in range(n_vertices, n_vertices + n)]

    omega_vertices = [list(range(n_vertices + n, n_vertices + n + n_omega))]

    if n_epsilon > 0: epsilon_vertices = [list(range(n_vertices + n + n_omega, 
                                                        n_vertices + n + n_omega + n_epsilon))]
    else: epsilon_vertices = []
    
    partition = inner_vertices + numbered_vertices + omega_vertices + epsilon_vertices

    return partition


def plot_graph(G, n_vertices, n, n_omega, n_epsilon):
    
    GG = Graph(G)  
    
    partition = get_partition(n_vertices, n, n_omega, n_epsilon)

    assert len(partition) >= n + 1

    color_dict = {}

    # inner vertices
    has_inner_vertices = (n_vertices > 0)
    if has_inner_vertices:
        color_dict.update({"gray":partition[0]}) 

    # numbered hairs
    cmap = plt.get_cmap('viridis')
    for i in range(1, n + has_inner_vertices): 
        assert len(partition[i]) == 1
        interpolation_value = (i - 1) / n
        color_dict.update({cmap(interpolation_value):partition[i]})

    # omega vertices
    color_dict.update({"red":partition[n + has_inner_vertices]}) 

    # epsilon vertices
    if n_epsilon > 0:
        color_dict.update({"magenta":partition[n + 1 + has_inner_vertices]}) 

    return GG.plot(vertex_colors=color_dict, vertex_labels=True)


def operate_on_EpsToOmega(G, n_vertices, n, n_omega, n_epsilon):
    """
    This is the "dual" operator to the delta_omega operator defined in https://arxiv.org/pdf/2407.16372 on page 4 
    which changes one omega vertex to an epsilon vertex.

    it operates on the graph G by:
    - for each epsilon vertex: replace it with an omega vertex and append the new graph to the image list
    - return the image list

    where replacing the epsilon with an omega vertex works as follows:
    - choose an epsilon vertex
    - note that the first epsilon vertex becomes an omega vertex since target.n_omega = domain.n_omega + 1
    - hence we simply swap the label of the chosen epsilon vertex with the label of the first epsilon vertex
    - the sign is given by the sign of the edge-permutation induced by the relabelling of the vertices
    - note that we do not need to additionaly consider the permutation of the omega-labels, since these are left in place by construction!

    EXAMPLE: "2x omaga & 5x epsilon attached to a single vertex"
    notation: "0(i)":
    - "i" denotes the "physical" vertex in the graph
    - "0" denotes the corresponding vertex-label
    1.  G1.vertices():  0(i), 1(o1), 2(o2), 3(e1), 4(e2), 5(e3), 6(e4), 7(e5)
    2.  eps_index = 5 -> eps_vertex e3
    3.  we now want to swap vertex e3 with e1: 
        relabeling_perm: [0, 1, 2, 5, 4, 3, 6, 7] 
        G2.vertices():  0(i), 1(o1), 2(o2), 3(e3), 4(e2), 5(e1), 6(e4), 7(e5)
    4.  Now, since each hair-vertex (o1, o2, e1, e2, e3, e4, e5) corresponds to exactly one edge (the hair),
        the induced edge permutation is simply a swap of (i-e3) with (i-e1)
        -> sign = -1
    """

    n_epsilon = len(G.vertices()) - n_vertices - n - n_omega
    assert n_epsilon >= 0

    G1 = copy(G)
    sgn = 1
    image = []

    # label all edges to determine sign later
    Shared.enumerate_edges(G1)

    for i in range(n_epsilon):
        
        omega_eps_cutoff = n_vertices + n + n_omega
        eps_index = omega_eps_cutoff + i
        
        # permute chosen epsilon vertex to the position where it becomes 
        # an omega vertex due to the new partition in self.domain, since:
        # self.target.n_omega = self.domain.n_omega + 1 
        # relabel: swap omega_eps_cutoff and eps_index
        # if i == 0: then omega_eps_cutoff == eps_index
        G2 = copy(G1)
        if i > 0:
            relabeling_perm = list(range(omega_eps_cutoff)) + [eps_index] + list(range(omega_eps_cutoff+1, eps_index)) + [omega_eps_cutoff] + list(range(eps_index+1, G1.order()))
            print(relabeling_perm)
            assert set(relabeling_perm) == set(G2.vertices())
            G2.relabel(relabeling_perm)

        # sign equals the perm_sign of the edge-permutation induced by the relabeling of the vertices
        sgn = Shared.shifted_edge_perm_sign2(G2)
        print("sgn: ", sgn)
        image.append((G2, sgn))

    return image


def plot_image_list(image_list, n_vertices, n, n_omega, n_epsilon):
    for i, image in enumerate(image_list):
        G, sgn = image
        #print(G.vertices())
        #print(G.edges())
        filename = "image_plot_%s.png" % str(i)
        filepath = "/mnt/c/Users/ps200/Desktop/Semesterarbeit/Plots/" + filename
        plot_graph(G, n_vertices, n, n_omega+1, n_epsilon-1).save(filepath)


def test_graph(G, n_vertices, n, n_omega, n_epsilon):

    filename = "_original_graph.png" 
    filepath = "/mnt/c/Users/ps200/Desktop/Semesterarbeit/Plots/" + filename
    plot_graph(G, n_vertices, n, n_omega, n_epsilon).save(filepath)
    
    image_list = operate_on_EpsToOmega(G, n_vertices, n, n_omega, n_epsilon)
    assert len(image_list) == n_epsilon
    plot_image_list(image_list, n_vertices, n, n_omega, n_epsilon)


# examples ---


# omega-epsilon double-leg
#test_graph(Graph({0:[1]}), n_vertices=0, n=0, n_omega=1, n_epsilon=1)

# 2x epsilon double-leg
#test_graph(Graph({0:[1]}), n_vertices=0, n=0, n_omega=0, n_epsilon=2)

# 1x omaga & 2x epsilon attached to a single vertex
# test_graph(Graph({0:[1, 2, 3]}), n_vertices=1, n=0, n_omega=1, n_epsilon=2)

# 2x omaga & 5x epsilon attached to a single vertex
test_graph(Graph({0:[1, 2, 3, 4, 5, 6, 7]}), n_vertices=1, n=0, n_omega=2, n_epsilon=5)
