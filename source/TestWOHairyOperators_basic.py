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


def operate_on_ContractEdges(G, n_vertices, n, n_omega, n_epsilon):

    n_old_eps = n_epsilon

    # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
    image = []
    for (i, e) in enumerate(G.edges(labels=False, sort=True)):

        print("contracting edge (u,v) =", e, "###########")
        
        (u, v) = e
        sgn = 1 

        # ensure u<v (this should be always true anyway actually)
        assert u < v

        # only edges connected to at least one internal vertex, and not connected to a numbered hair-vertex can be contracted

        # both u and v are not internal vertices
        if u >= n_vertices: continue

        # u or v are numbered vertices        
        if u >= n_vertices and u < n_vertices + n: continue
        if v >= n_vertices and v < n_vertices + n: continue

        # label all edges to determine sign later
        G1 = copy(G)
        Shared.enumerate_edges(G1)
        print("intial state: ")
        print("vertices: ", G1.vertices())
        print("edges: ", G1.edges())


        # contracting the edge (u,v)
        # we always delete the lower index vertex. This ensures that the extra vertices are never deleted

        # if both u and v are internal vertices
        if False and v < n_vertices:
            assert u < n_vertices

            # contract edge by merging vertices
            G1.merge_vertices([v, u])

            # if we have too few edges, some double edges have been created => zero
            if G.size() - G1.size() > 1: continue

            # relabel vertices back to 0,1,2,...,k
            G1.relabel(range(len(G1.vertices())), inplace=True)
            
            # find edge permutation sign
            sgn *= Shared.shifted_edge_perm_sign2(G1)
            # print("sgn3_",sgn)
            print("adding to image ---")
            print("vertices: ", G1.vertices())
            print("edges: ", G1.edges())
            image.append((G1, sgn))
            # image.append((Graph(G1.graph6_string()), sgn))
            # print("hmm0:", G.graph6_string(), G1.graph6_string())


        # if u is an internal vertex and v is an epsilon-vertex
        elif False and u < n_vertices \
            and v >= n_vertices + n + n_omega:
            
            G1.delete_vertex(v)

            # split u into separate epsilon-vertices
            for w in G1.neighbors(u):

                print("picking neighbour w =", w, " ---")
                
                u_w_label = G1.edge_label(u, w)
                G1.delete_edge(u, w)

                new_eps_label = max(G1.vertices()) + 1
                assert not new_eps_label in G1.vertices()

                G1.add_vertex(new_eps_label)
                G1.add_edge(w, new_eps_label, u_w_label)

            G1.delete_vertex(u) 

            # relabel vertices back to 0,1,2,...,k
            G1.relabel(range(len(G1.vertices())), inplace=True)

            sgn *= Shared.shifted_edge_perm_sign2(G1)
            print("adding to image ---")
            print("sign:", sgn)
            print("vertices: ", G1.vertices())
            print("edges: ", G1.edges())
            image.append((G1, sgn))


        

        # if u is an internal vertex and v is a omega-vertex
        # the second vertex is now an omega-vertex, so we need to merge the vertex with the eps vertex
        # after reconnecting one of the edges to omega
        # we assume that u != eps, because eps-omega-edges cannot be contracted
        if True: continue
        elif u < n_vertices \
            and v >= n_vertices + n \
            and v < n_vertices + n + n_omega:
            
            pass
            # TODO

    return image

"""
def operate_on_ContractEdges(G, n_vertices, n, n_omega, n_epsilon):

    n_old_eps = n_epsilon

    # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
    image = []
    for (i, e) in enumerate(G.edges(labels=False, sort=True)):

        print("contracting edge (u,v) =", e)
        
        (u, v) = e
        sgn = (-1)**i

        # ensure u<v (this should be always true anyway actually)
        assert u < v

        # only edges connected to at least one internal vertex, and not connected to a numbered hair-vertex can be contracted

        # both u and v are not internal vertices
        if u >= n_vertices: continue

        # u or v are numbered vertices        
        if u >= n_vertices and u < n_vertices + n: continue
        if v >= n_vertices and v < n_vertices + n: continue
    
        # print("sgn0",sgn)
        previous_size = G.size()

        # print("sgn1",sgn)
        G1 = copy(G)
        # label all edges to determine sign later
        Shared.enumerate_edges(G1)
        print("intial state: ###########")
        print("vertices: ", G1.vertices())
        print("edges: ", G1.edges())


        # contracting the edge (u,v)
        # we always delete the lower index vertex. This ensures that the extra vertices are never deleted

        # if both u and v are internal vertices
        if v < n_vertices:
            assert u < n_vertices

            # contract edge by merging vertices
            G1.merge_vertices([v, u])

            # if we have too few edges some double edges have been created => zero
            if (previous_size - G1.size()) > 1: continue

            # relabel vertices
            G1.relabel(range(len(G1.vertices())), inplace=True)
            
            # find edge permutation sign
            sgn *= Shared.shifted_edge_perm_sign2(G1)
            # print("sgn3_",sgn)
            print("adding to image ---")
            print("vertices: ", G1.vertices())
            print("edges: ", G1.edges())
            image.append((G1, sgn))
            # image.append((Graph(G1.graph6_string()), sgn))
            # print("hmm0:", G.graph6_string(), G1.graph6_string())


        # if u is an internal vertex and v is a omega-vertex
        # the second vertex is now an omega-vertex, so we need to merge the vertex with the eps vertex
        # after reconnecting one of the edges to omega
        # we assume that u != eps, because eps-omega-edges cannot be contracted
        elif u < n_vertices \
            and v >= n_vertices + n \
            and v < n_vertices + n + n_omega:
            
            u_v_label = G1.edge_label(u, v)
            G1.delete_edge(u, v)
            print("state after deleting u_v:")
            print("vertices: ", G1.vertices())
            print("edges: ", G1.edges())

            # loop over neighbors w of u to be connected to omega
            for w in G1.neighbors(u):

                print("picking neighbour w =", w, " ---")

                G2 = copy(G1)
                sgn2 = sgn

                # u_w edge -> v_w edge
                u_w_label = G2.edge_label(u, w)
                G2.delete_edge(u, w)

                G2.add_edge(v, w, min(u_v_label, u_w_label))

                print("state after u_w -> v_w:")
                print("vertices: ", G2.vertices())
                print("edges: ", G2.edges())

                # we now need to split u into separate epsilon-vertices
                # for each remaining adjecent vertex s of u
                n_new_eps = len(G2.neighbors(u))
                for (j, s) in enumerate(G2.neighbors(u)): 
                    
                    u_s_label = G2.edge_label(u, s)
                    G2.delete_edge(u, s)

                    # add new epsilon-vertex
                    # since we add the new vertex at the end, it automatcally corresponds to epsilon
                    assert set(G2.vertices()) == set(range(len(G2.vertices())))
                    new_eps_label = len(G2.vertices()) 

                    old_size = G2.order()
                    #print(G2.vertices())
                    G2.add_vertex(new_eps_label)
                    assert G2.order() == old_size + 1, (G2.order(), old_size)

                    if j == 0: label = min(max(u_v_label, u_w_label), u_s_label)
                    else: label = len(G2.edges()) + (j-1)
                    G2.add_edge(new_eps_label, s, label)

                G2.delete_vertex(u)

                print("state after turning u into epsilon(s):")
                print("vertices: ", G2.vertices())
                print("edges: ", G2.edges())

                # in case we have too few edges some double edges have been created => zero
                if (previous_size - G2.size()) > 1: continue

                G2.relabel(range(len(G2.vertices())), inplace=True)
                # find edge permutation sign
                sgn2 *= Shared.shifted_edge_perm_sign2(G2)
                
                # sanity checks
                #n_eps_in_target = self.target.get_n_epsilon_from_graph(G2)
                #assert n_eps_in_target == n_old_eps + n_new_eps
                #assert G2.order() == self.target.n_vertices + self.target.n + self.target.n_omega + n_eps_in_target
                print("adding to image ---")
                print("sign:", sgn2)
                print("vertices: ", G1.vertices())
                print("edges: ", G1.edges())
                image.append((G2, sgn2))

    return image
"""

def plot_image_list(image_list, n_vertices, n, n_omega, n_epsilon, operator):
    for i, image in enumerate(image_list):
        G, sgn = image
        #print(G.vertices())
        #print(G.edges())
        filename = "image_plot_%s.png" % str(i)
        filepath = "/mnt/c/Users/ps200/Desktop/Semesterarbeit/Plots/" + filename

        if operator == operate_on_EpsToOmega:
            plot_graph(G, n_vertices, n, n_omega+1, n_epsilon-1).save(filepath)

        elif operator == operate_on_ContractEdges:
            n_vertices_new = n_vertices - 1
            n_epsilon_new = G.order() - n_vertices_new - n - n_omega
            assert n_vertices_new >= 0 and n_epsilon_new >= 0
            plot_graph(G, n_vertices_new, n, n_omega, n_epsilon_new).save(filepath)


def test_graph(G, n_vertices, n, n_omega, n_epsilon, operator):

    print("testing graph: ", G, operator)
    assert G.order() == n_vertices + n + n_omega + n_epsilon

    filename = "_original_graph.png" 
    filepath = "/mnt/c/Users/ps200/Desktop/Semesterarbeit/Plots/" + filename
    plot_graph(G, n_vertices, n, n_omega, n_epsilon).save(filepath)
    
    image_list = operator(G, n_vertices, n, n_omega, n_epsilon)
    print("image_list: ", image_list)
    plot_image_list(image_list, n_vertices, n, n_omega, n_epsilon, operator)

    # check number of image-graphs
    if operator == operate_on_EpsToOmega:
        assert len(image_list) == n_epsilon



# Testing ContractEdges ---

# 1x omaga & 2x epsilon attached to a single vertex
# test_graph(Graph({0:[1, 2, 3]}), n_vertices=1, n=0, n_omega=1, n_epsilon=2, operator=operate_on_ContractEdges)

# 4x eps attached to two vertices connected by an internal edge, which is being contracted -> 4x eps attached to vertex
#test_graph(Graph({0:[1, 2, 3], 1:[4,5]}), n_vertices=2, n=0, n_omega=0, n_epsilon=4, operator=operate_on_ContractEdges)


test_graph(Graph({0:[1, 2, 3], 1:[2,4]}), n_vertices=3, n=0, n_omega=0, n_epsilon=2, operator=operate_on_ContractEdges)




# Testing EpsToOmega ---

# omega-epsilon double-leg
#test_graph(Graph({0:[1]}), n_vertices=0, n=0, n_omega=1, n_epsilon=1, operator=operate_on_EpsToOmega)

# 2x epsilon double-leg
#test_graph(Graph({0:[1]}), n_vertices=0, n=0, n_omega=0, n_epsilon=2, operator=operate_on_EpsToOmega)

# 1x omaga & 2x epsilon attached to a single vertex
# test_graph(Graph({0:[1, 2, 3]}), n_vertices=1, n=0, n_omega=1, n_epsilon=2, operator=operate_on_EpsToOmega)

# 2x omaga & 5x epsilon attached to a single vertex
#test_graph(Graph({0:[1, 2, 3, 4, 5, 6, 7]}), n_vertices=1, n=0, n_omega=2, n_epsilon=5, operator=operate_on_EpsToOmega)

