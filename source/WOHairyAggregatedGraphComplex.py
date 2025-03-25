"""
Graph complexes based on simple graphs with numbered hairs and hairs of two (omega- and epsilon-)decorations
as in the graph complex computing weight 11 cohomology. The omega decorations are considered odd.
For more detailed Information see the following references:
- Weight 11 compactly supported cohomology of moduli spaces of curves: https://arxiv.org/abs/2302.04204
- Weight 11 Compactly Supported Cohomology of Moduli Spaces of Curves in excess four: https://arxiv.org/abs/2407.16372

This is implemented in the following files:
    _WOHairyComponentsGraphComplex.py 
    _WOHairyAggregatedGraphComplex.py <--- (this file)
    _WOHairyFinalGraphComplex.py

This specific file _WOHairyAggregatedGraphComplex.py recursively joins together connected graphs 
from _WOHairyComponentsGraphComplex.py to get basis-graphs with any number of connected components.
"""



from sage.all import *
import Shared
import Parameters
import WOHairyComponentGraphComplex
import OrdinaryGraphComplex
import itertools
import math

graph_type = "wohairy_aggregated"

# ------- Graph Vector Space --------
class WOHairyAggregatedGVS(WOHairyComponentGraphComplex.WOHairyComponentGVS):

    def __init__(self, n_components, n_vertices, genus, n, n_omega, n_epsilon, n_double_legs, do_print=False):

        self.n_components = n_components
        self.n = n
        self.n_vertices = n_vertices
        self.genus = genus
        self.n_omega = n_omega
        self.n_epsilon = n_epsilon
        self.n_double_legs = n_double_legs
        self.n_hairs = self.n + self.n_omega + self.n_epsilon
        self.total_num_vertices = self.n_vertices + self.n_hairs
        self.excess = 3*(self.genus - 1) + 2*self.n - 2*self.n_omega

        # computation of the inner edges by genus
        # g = (E-V+num_components) + 1 of the contracted graph:
        # g = (n_edges + n_hairs - n_double_legs) - (n_vertices + n + 1) + 1 + 1 
        #   = n_edges + n_hairs - n_vertices - n - n_double_legs + num_components
        self.n_edges = genus + n_vertices + n + n_double_legs - self.n_hairs - 1

        # loop-order of the inner graph
        self.n_loops = self.n_edges - n_vertices + (self.n_components - self.n_double_legs)

        # corresponding ogvs for perm_sign method
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(self.n_vertices + self.n_hairs, self.n_loops, even_edges=False)

        if do_print:
            print("intializing WOHairyAggregatedGVS with params: -------")
            print("n_components: ", n_components)
            print("n_vertices: ", n_vertices)
            print("genus: ", genus)
            print("n: ", n)
            print("n_omega: ", n_omega)
            print("n_epsilon: ", n_epsilon)
            print("n_double_legs: ", n_double_legs)
            print("n_edges: ", self.n_edges)
            print("n_loops:", self.n_loops)
            print("is_valid: ", str(self.is_valid()))


    def cohom_degree(self):
        # deg = 22 + (all_edges - n_omega - n)
        # where all_edges = inner_edges + n + n_eps + n_omega - n_double_legs
        # -> deg = 22 + inner_edges + n_eps - n_double_legs
        return 22 + self.n_edges + self.n_epsilon - self.n_double_legs
    

    def __hash__(self):
        return hash("wo_agg_gra%d_%d_%d_%d_%d_%d_%d" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "wo_agg_gra%d_%d_%d_%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_basis_file_path(self):
        s = "wo_agg_gra%d_%d_%d_%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('components', self.n_components), ('vertices', self.n_vertices), ('genus', self.genus), ('numbered_hairs', self.n), 
                                   ('omegas', self.n_omega), ('epsilons', self.n_epsilon), ('double_legs', self.n_double_legs)])


    def is_valid(self):
        # possible values
        l = self.n_components >= 1 \
        and self.n_vertices >= 0 \
        and self.genus >= 1 \
        and self.n >= 0 \
        and self.n_omega >= 0 \
        and self.n_epsilon >= 0 \
        and self.n_double_legs >= 0 \
        and self.n_hairs >= 0 \
        and self.n_edges >= 0 \
        and self.n_loops >= 0 \
        and self.excess >= 0
        #print("is_valid() after value-checking:", l)
        # each n_double_leg corresponds to a connected component 
        l = l and self.n_double_legs <= self.n_components
        # number of components which are not double-legs: self.n_components - self.n_double_legs
        l = l and self.n_vertices >= self.n_components - self.n_double_legs
        # If all components are double-legs:
        if self.n_components == self.n_double_legs: 
            l = l and self.n_vertices == 0
            l = l and self.n_hairs == 2*self.n_double_legs
            l = l and self.n_loops == 0
            l = l and self.n_edges == 0
        
        # each connected component must at least contain one omega or epsilon
        l = l and self.n_omega + self.n_epsilon >= self.n_components
        # At least trivalent internal vertices. (each double_leg removes two hair-decorations from the equation)
        l = l and (3 * self.n_vertices <= 2 * self.n_edges + (self.n_hairs - 2*self.n_double_legs)) 
        # At most a full graph.
        l = l and self.n_edges <= (self.n_vertices) * (self.n_vertices - 1) / 2
        # each double-leg contains exactly two hair lables (which are not inner vertices)
        l = l and self.n_double_legs <= (self.total_num_vertices - self.n_vertices) / 2
        # each other connected component contains at least one vertex
        l = l and self.n_components - self.n_double_legs <= self.total_num_vertices - 2*self.n_double_legs
        # due to symmetry reasons, each vertex can have maximally one epsilon-hair attached to it otherwise the graph is equal to zero
        l = l and self.n_epsilon <= 2*self.n_double_legs + self.n_vertices
        # cannot have more eps than excess
        l = l and self.n_epsilon <= self.excess
        
        # test for single-vertex trees with omega and numbered hairs
        if self.n_vertices == 1 and self.n_components == 1 and self.n_epsilon == 0:
            l = l and self.n_double_legs == 0
            l = l and (self.excess - self.n_omega + 3) % 2 == 0
            l = l and self.excess == self.n_omega - 3 + 2*self.n
            l = l and self.n_omega + self.n >= 3
        

        # lemma: for excess <= 4, the loop order is 0
        if self.excess <= 4: l = l and self.n_loops == 0
        
        return l


    

    def get_generating_graphs(self):
        #print("Aggregated: get_generating_graphs")
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            print("self is not valid")
            return []
    

        # recursive aggregation by n_connected components
        assert isinstance(self.n_components, int)
        assert (self.n_components >= 1)

        # if we only have one connected component the definition is equal to the definition of ModifedCHairyGraphComplex
        if self.n_components == 1:

            V = WOHairyComponentGraphComplex.WOHairyComponentGVS(n_vertices=self.n_vertices,
                                                              n_loops=self.n_loops,
                                                              n=self.n,
                                                              n_omega=self.n_omega,
                                                              n_epsilon=self.n_epsilon)

            V.build_basis(ignore_existing_files=False)

            #print("prebuilt basis len:", len([G for G in V.get_basis()]))

            for G in V.get_basis():

                assert len(G.vertices()) == self.total_num_vertices
                assert len(G.edges()) == self.n_edges + self.n_hairs - self.n_double_legs

                yield G


        # else we need to iterate through different possible configurations of joining 
        # a graph complex with n_components-1 connected components with one with only a single connected component.
        else: 
            #print("recursion-step for n_components=" + str(self.n_components))

            configurations = self.get_configurations()
            #print(configurations)
            #print("number of combinatorically possible configurations: " + str(len(configurations)))
            #print(configurations)
            for configuration in configurations:
                #print("configuration: " + str(configuration))
                (n_vert_1, n_vert_2, genus_1, genus_2, n_1, n_2, n_omega_1, n_omega_2, n_eps_1, n_eps_2, n_DL_1, n_DL_2) = configuration

                n_components_1 = 1
                n_components_2 = self.n_components - 1
                
                V_1 = WOHairyAggregatedGVS(n_components=n_components_1,
                                         n_vertices=n_vert_1,
                                         genus=genus_1,
                                         n=n_1,
                                         n_omega=n_omega_1,
                                         n_epsilon=n_eps_1,
                                         n_double_legs=n_DL_1)
                
                if not V_1.is_valid(): continue

                V_2 = WOHairyAggregatedGVS(n_components=n_components_2,
                                         n_vertices=n_vert_2,
                                         genus=genus_2,
                                         n=n_2,
                                         n_omega=n_omega_2,
                                         n_epsilon=n_eps_2,
                                         n_double_legs=n_DL_2)
                
                if not V_2.is_valid(): continue

                # check if excess adds up
                if self.excess != V_1.excess + V_2.excess: continue

                V_1.build_basis(ignore_existing_files=False)
                #print("built basis 1")
                V_2.build_basis(ignore_existing_files=False)
                #print("built basis 2")
                #print("dimension V_1:", V_1.get_dimension())
                #print("dimension V_2:", V_2.get_dimension())

                if V_1.get_dimension() > 0 and V_2.get_dimension()>0 and False:
                    print("actually possible configuration", configuration)

                    print("n_components_1", n_components_1)
                    print("n_vert_1", n_vert_1)
                    print("genus_1", genus_1)
                    print("n_1", n_1)
                    print("n_omega_1", n_omega_1)
                    print("n_eps_1", n_eps_1)
                    print("n_DL_1", n_DL_1)
                
                    print("n_components_2", n_components_2)
                    print("n_vert_2", n_vert_2)
                    print("genus_2", genus_2)
                    print("n_2", n_2)
                    print("n_omega_2", n_omega_2)
                    print("n_eps_2", n_eps_2)
                    print("n_DL_2", n_DL_2)
                
                
                # join graphs
                #print("joining possible graph combinations")
                for G_single_comp in V_1.get_basis():
                    for G2_mult_comp in V_2.get_basis():
                        #print("joining graphs")

                        # join graphs:
                        # G now has vertices (0,i) and (1,j)
                        G = G_single_comp.disjoint_union(G2_mult_comp)

                        assert len(G.vertices()) == self.total_num_vertices
                        assert len(G.edges()) == self.n_edges + self.n_hairs - self.n_double_legs

                        # relabel vertices to 0,1,2,...
                        G.relabel({old: new for new, old in enumerate(G.vertices())}, inplace=True)
                
                        # permute vertices accordingly to how the partition is made
                        old_order = G.vertices()
                        new_orders = reorder_vertices(old_order, n_vert_1, n_vert_2, n_1, n_2, n_omega_1, n_omega_2, n_eps_1, n_eps_2)

                        for new_order in new_orders:
                            mapping = {new:old for old, new in zip(old_order, new_order)}
                            
                            G_new = G.copy()
                            G_new.relabel(mapping)

                            yield G_new

                            """
                            print("G.vertices(): " + str(G.vertices()))
                            print("G.edges(): " + str(G.edges()))
                            print("H.vertices(): " + str(H.vertices()))
                            print("H.edges(): " + str(H.edges()))
                            """
        
    


    def get_configurations(self):
        
        configurations = []
        
        for n_vert_1 in range(0, self.n_vertices + 1):
            n_vert_2 = self.n_vertices-n_vert_1
            #print("n_vert_1: " + str(n_vert_1))

            # ensure that genus_1 and genus_2 are both >= 1
            for genus_1 in range(1, self.genus+1):
                # g = (E-V+1) + 1 of the contracted graph:
                #   = (n_edges + n_hairs - n_double_legs) - (n_vertices + n + 1) + 1 + 1 
                #   = n_edges + n_hairs - n_vertices - n_double_legs - n + 1 
                # g_tot = g_1 + g_2 - 1
                genus_2 = self.genus - genus_1 + 1
                assert genus_1 >= 1 and genus_2 >= 1
                #print("genus_1: " + str(genus_1))

                for n_1 in range(self.n + 1):
                    n_2 = self.n-n_1
                    #print("n_1: " + str(n_1))

                    for n_omega_1 in range(self.n_omega + 1):
                        n_omega_2 = self.n_omega-n_omega_1
                        #print("n_omega_1: " + str(n_omega_1))

                        for n_epsilon_1 in range(self.n_epsilon + 1):
                            n_epsilon_2 = self.n_epsilon-n_epsilon_1
                            #print("n_epsilon_1: " + str(n_epsilon_1))

                            for n_double_legs_1 in range(0, self.n_double_legs+1):
                                n_double_legs_2 = self.n_double_legs - n_double_legs_1
                            
                                if (n_omega_1 >= 1 or n_epsilon_1 >= 1) and (n_omega_2 >= 1 or n_epsilon_2 >= 1):

                                    configurations.append((n_vert_1, n_vert_2,
                                                        genus_1, genus_2,
                                                        n_1, n_2,
                                                        n_omega_1, n_omega_2,
                                                        n_epsilon_1, n_epsilon_2,
                                                        n_double_legs_1, n_double_legs_2))
        
        return configurations 

        


# helper-functions ---
        

def list_elems_from_range(input_list, num, dist):
    return [input_list[i] for i in list(range(num, num + dist))]

def list_elems_from_double_range(input_list, num1, dist1, num2, dist2):
    """
    input: some list 
    output: two disjoint slices of the list (starting at num_i and of length dist_i)
            merged together to a single new list
    """
    assert max(num1 + dist1, num2 + dist2) <= len(input_list)
    assert num1 + dist1 <= num2 or num2 + dist2 <= num1 # disjointness of the slice

    return list_elems_from_range(input_list, num1, dist1) + list_elems_from_range(input_list, num2, dist2)


def get_cross_permutations(list_1, list_2):
    """
    Since the individual components already contain all of the internal permutations of the numbered hairs, 
    we only need to consider order-preserving permutations which swap numbered hairs between the two components.
    Let n_1, n_2 be the respective number of numbered hairs
    Then within the generation of the individual components we have already considered 
    the n_i! permutations of the elements in list_i for i=1,2.
    Hence in order to consider all (n_1+n_2)! permutations it is only left to check 
    additional (n_1+n_2)! / (n_1!*n_2!) permutations for any combination of such components.

    More operationally we compute the following:
    Generate all permutations of S that satisfy the constraints:
      - S is partitioned into disjoint subsets A and B.
      - Each element in A is either fixed or swapped with an element in B (and vice versa).
      - If two elements i, j in A (i < j) are swapped, their images in B preserve order (σ(i) < σ(j)).
        Likewise, if two in B are swapped, their images in A preserve order.
    Yields each valid permutation as a tuple of length n (in one-line notation).
    """
    
    # assert that list_1 and list_2 are disjoint with n combined elements
    assert len(list_1) + len(list_2) == len(set(list_1 + list_2))

    permutations = []

    # Iterate over possible swap counts
    for k in range(0, min(len(list_1), len(list_2)) + 1):
        # Choose k elements from A and k from B to swap
        for combA in itertools.combinations(list_1, k):
            #print("combA", combA)
            for combB in itertools.combinations(list_2, k):
                #print("combB", combB)
                # Pair the chosen elements in sorted order (combinations are already sorted)
                perm = list_1 + list_2  # start with identity permutation
                for a_elem, b_elem in zip(combA, combB):
                    # Swap a_elem and b_elem in the permutation
                    #print("before swap", perm)
                    i = perm.index(a_elem)
                    #print(i)
                    j = perm.index(b_elem)
                    #print(j)
                    
                    # Swap in place
                    perm[i], perm[j] = perm[j], perm[i]
                    #print("after swap", perm)
                    
                permutations.append(perm)
    
    
    n_1 = len(list_1)
    n_2 = len(list_2)
    assert len(permutations) == math.factorial(n_1 + n_2) / (math.factorial(n_1) * math.factorial(n_2))

    return permutations



def reorder_vertices(old_order, n_vert_1, n_vert_2, n_1, n_2, n_omega_1, n_omega_2, n_epsilon_1, n_epsilon_2):
    assert len(old_order) == n_vert_1 + n_vert_2 + n_1 + n_2 + n_omega_1 + n_omega_2 + n_epsilon_1 + n_epsilon_2

    new_orders = []
    #print("old_order: " + str(old_order))
    divider = n_vert_1 + n_1 + n_omega_1 + n_epsilon_1

    inner_vertices = list_elems_from_double_range(old_order, 0, n_vert_1, divider, n_vert_2)
    #print("inner_vertices: " + str(inner_vertices))
    
    numbered_vertices_1 = list_elems_from_range(old_order, n_vert_1, n_1)
    numbered_vertices_2 = list_elems_from_range(old_order, divider+n_vert_2, n_2)
    #print("numbered_vertices: " + str(numbered_vertices))

    omega_vertices = list_elems_from_double_range(old_order, n_vert_1+n_1, n_omega_1, divider+n_vert_2+n_2, n_omega_2)
    #print("omega_vertices: " + str(omega_vertices))

    epsilon_vertices = list_elems_from_double_range(old_order, n_vert_1+n_1+n_omega_1, n_epsilon_1, divider+n_vert_2+n_2+n_omega_2, n_epsilon_2)
    #print("epsilon_vertices: " + str(epsilon_vertices))

    #print("list1", numbered_vertices_1)
    #print("list2", numbered_vertices_2)
    #for permuted_numbered_verticies in itertools.permutations(numbered_vertices_1 + numbered_vertices_2):
    #print("un-permuted:", numbered_vertices_1 + numbered_vertices_2)
    for permuted_numbered_verticies in get_cross_permutations(numbered_vertices_1, numbered_vertices_2):
        
        #print("permuted_numbered_verticies: " + str(permuted_numbered_verticies))
        new_order = inner_vertices + list(permuted_numbered_verticies) + omega_vertices + epsilon_vertices

        #print("new_order: " + str(new_order))
        assert len(new_order) == len(old_order)
        assert set(new_order) == set(range(len(new_order)))

        new_orders.append(new_order)

    return new_orders

