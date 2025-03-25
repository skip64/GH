"""
Graph complexes based on simple graphs with numbered hairs and hairs of two (omega- and epsilon-)decorations
as in the graph complex computing weight 11 cohomology. The omega decorations are considered odd.
For more detailed Information see the following references:
- Weight 11 compactly supported cohomology of moduli spaces of curves: https://arxiv.org/abs/2302.04204
- Weight 11 Compactly Supported Cohomology of Moduli Spaces of Curves in excess four: https://arxiv.org/abs/2407.16372

This is implemented in the following files:
    _WOHairyComponentsGraphComplex.py <--- (this file)
    _WOHairyAggregatedGraphComplex.py
    _WOHairyFinalGraphComplex.py

This specific file _WOHairyComponentGraphComplex.py generates the connected compontents of such graphs.
This is achieved as follows:
1. generate Hairy graphs with numbered hairs using CHairyGraphComplex.py
2. relable some of the hairs with omega or epsilon
"""



from sage.all import *
import Shared
import Parameters
import CHairyGraphComplex
import itertools
import OrdinaryGraphComplex

graph_type = "wohairy_component"


# ------- Graph Vector Space --------
class WOHairyComponentGVS(CHairyGraphComplex.CHairyGraphVS):

    def __init__(self, n_vertices, n_loops, n, n_omega, n_epsilon):

        self.n = n
        self.n_vertices = n_vertices
        self.n_loops = n_loops # inner loops! (without merging epsilons and omegas)
        self.n_omega = n_omega
        self.n_epsilon = n_epsilon
        self.even_edges = False
        self.n_hairs = self.n + self.n_omega + self.n_epsilon

        # computation of inner-edges
        if self.n_vertices > 0: self.n_edges = self.n_loops + self.n_vertices - 1 # case of no double-leg
        else: self.n_edges = 0 # case of double-leg

        # corresponding ogvs for perm_sign method
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(self.n_vertices + self.n_hairs, self.n_loops, even_edges=False)



    def __hash__(self):
        return hash("wo_comp_gra%d_%d_%d_%d_%d" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "wo_comp_gra%d_%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_basis_file_path(self):
        s = "wo_comp_gra%d_%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops), ('numbered_hairs', self.n), 
                                   ('omegas', self.n_omega), ('epsilons', self.n_epsilon)])


    def get_partition(self):

        if self.n_vertices > 0: inner_vertices = [list(range(0, self.n_vertices))]
        else: inner_vertices = []

        numbered_vertices = [[j] for j in range(self.n_vertices, self.n_vertices + self.n)]

        omega_vertices = [list(range(self.n_vertices + self.n, self.n_vertices + self.n + self.n_omega))]

        if self.n_epsilon > 0: epsilon_vertices = [list(range(self.n_vertices + self.n + self.n_omega, 
                                                            self.n_vertices + self.n + self.n_omega + self.n_epsilon))]
        else: epsilon_vertices = []
        
        partition = inner_vertices + numbered_vertices + omega_vertices + epsilon_vertices

        return partition



    def is_valid(self):

        # allowed values
        l = self.n_vertices >= 0 and self.n_loops >= 0 and self.n_edges >= 0 \
        and self.n >= 0 and self.n_omega >= 0 and self.n_epsilon >= 0
        # At least trivalent internal vertices.
        l = l and (3 * self.n_vertices <= 2 * self.n_edges + self.n_hairs)
        # At most a full graph.
        l = l and self.n_edges <= (self.n_vertices) * (self.n_vertices - 1) / 2
        # connected graph must contain at least one omega or epsilon
        l = l and self.n_omega >= 1 or self.n_epsilon >= 1
        # if double-leg
        if self.n_vertices == 0: l = l and self.n_hairs == 2

        return l
    

    def get_generating_graphs(self):

        if not self.is_valid():
            print("self is not valid")
            return []
    
        # produce all hairy graphs
        hairy_graphs = self.get_hairy_graphs(self.n_vertices, self.n_loops, self.n_hairs, include_novertgraph=True)

        # produce all neccesary permutations of the hairs
        # for more information look at the function "multiset_permutations" below
        all_perm = multiset_permutations(n_vertices=self.n_vertices, n=self.n, n_omega=self.n_omega, n_epsilon=self.n_epsilon)

        return (G.relabel(p, inplace=False) for G in hairy_graphs for p in all_perm)



    def perm_sign(self, G, p):
        # The sign is the same as the corresponding sign in the
        # ordinary graph complex, apart from an extra contribution from the omega-hair-vertices.

        sgn = self.ogvs.perm_sign(G, p)

        if self.n_omega > 0:
            omega_hairs = p[self.n_vertices + self.n : self.n_vertices + self.n + self.n_omega]
            sgn_omega_perm = Shared.Perm.shifted(omega_hairs).signature()
            sgn *= sgn_omega_perm

        return sgn
    

# helper-functions ---

def multiset_permutations(n_vertices, n, n_omega, n_epsilon):
    """
    this is a helper-function for WOHairyComponentGVS.get_generating_graphs()
    instead of considering all permutations of the hairs of the graph, we can make the following restriction:
    - omega-hairs are identical and hence need not be permuted with each other
    - the same holds for the epsilon-hairs

    Using this we need only to consider (n + n_omega + n_epsilon)! / (n_omega! * n_epsilon!)
    instead of (n + n_omega + n_epsilon)! permutations
    
    Example: ---
    - input: [n_vertices=0, n=1, n_omega=2, n_epsilon=0]
    - possible permutations of the hairs: [1, omega, omega], [omega, 1, omega], [omega, omega, 1]
    - intuetively the output should be: [0, 1, 2], [1, 0, 2], [2, 0, 1]
    - BUT: the partition is always given by [[0], [1,2]]
      hence the new "1"-vertex needs to get relabled to 0
    - because of this the output consits of the inverted permutations of the intuitive output:
        [0, 1, 2], [1, 0, 2], [1, 2, 0]

    """

    permutations = []

    L = n + n_omega + n_epsilon
    
    # 1) Choose positions for the n_omega identical omega's
    for omega_positions in itertools.combinations(range(L), n_omega):
        set_omega = set(omega_positions)
        
        # 2) Choose positions (from the leftover) for the n_epsilon identical epsilon's
        leftover_after_omega = [i for i in range(L) if i not in set_omega]
        for epsilon_positions in itertools.combinations(leftover_after_omega, n_epsilon):
            set_epsilon = set(epsilon_positions)
            
            # 3) The remaining L - n_omega - n_epsilon positions are for the n distinct items.
            leftover_for_distinct = [i for i in range(L)
                                     if i not in set_omega and i not in set_epsilon]
            
            # Permute the n distinct items in all possible ways
            for perm in itertools.permutations(leftover_for_distinct):
                
                # compute the full permutation
                full_perm = list(perm) + list(omega_positions) + list(epsilon_positions) 
                assert set(full_perm) == set(list(range(L)))

                # invert the permutation
                full_perm_inv = [None]*L
                for i in range(L):
                    full_perm_inv[full_perm[i]] = i
                
                # add the identity-permutaion for the inner vertices and
                # translate the permutation by the number of inner vertices 
                full_perm_inv_translated = list(range(0, n_vertices)) + [n_vertices + i for i in full_perm_inv] 

                permutations.append(full_perm_inv_translated)

    # assert that we have generated the correct number of permutations
    assert len(permutations) == math.factorial(n + n_omega + n_epsilon) / (math.factorial(n_omega) * math.factorial(n_epsilon))

    return permutations


    
