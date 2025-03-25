"""
Graph complexes based on simple graphs with numbered hairs and hairs of two (omega- and epsilon-)decorations
as in the graph complex computing weight 11 cohomology. The omega decorations are considered odd.
For more detailed Information see the following references:
- Weight 11 compactly supported cohomology of moduli spaces of curves: https://arxiv.org/abs/2302.04204
- Weight 11 Compactly Supported Cohomology of Moduli Spaces of Curves in excess four: https://arxiv.org/abs/2407.16372

This is implemented in the following files:
    _WOHairyComponentsGraphComplex.py 
    _WOHairyAggregatedGraphComplex.py
    _WOHairyFinalGraphComplex.py <--- (this file)

This specific file _WOHairyFinalGraphComplex.py selects, for some given cohomoligcal degree,
the relevant graphs from _WOHairyAggregatedGraphComplex.py and puts them in a single file.
"""



from sage.all import *
import Shared
import Parameters
import WOHairyAggregatedGraphComplex
import GraphVectorSpace
import matplotlib.pyplot as plt


graph_type = "wohairy_final"


# ------- Graph Vector Space --------
class WOHairyFinalGVS(GraphVectorSpace.GraphVectorSpace):

    def __init__(self, genus, n, n_omega, degree):

        # -> genus & degreee bestimmt anzahl vertices !
        # epsilons zusammenkleben sorgt für einheitliche Partition

        #print("initializing WOHairyFinalGVS ---")

        self.genus = genus
        self.n = n
        self.n_omega = n_omega
        self.degree = degree   
        
        self.n_vertices = degree - 22 + n_omega - genus + 1
        self.excess = 3*genus + 2*n - 25 


    def __hash__(self):
        return hash("wo_fin_gra%d_%d_%d_%d" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "wo_fin_gra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_basis_file_path(self):
        s = "wo_fin_gra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('genus', self.genus), ('n', self.n), 
                                   ('omegas', self.n_omega), ('degree', self.degree)])


    # since we don't a priori know the number of epsilons, we also need to pass the graph G
    def get_partition(self, G):
        
        n_epsilon = len(G.vertices()) - self.n_vertices - self.n - self.n_omega
        assert n_epsilon >= 0

        if self.n_vertices > 0: inner_vertices = [list(range(0, self.n_vertices))]
        else: inner_vertices = [[]]

        numbered_vertices = [[j] for j in range(self.n_vertices, self.n_vertices + self.n)]

        omega_vertices = [list(range(self.n_vertices + self.n, self.n_vertices + self.n + self.n_omega))]

        if n_epsilon > 0: epsilon_vertices = [list(range(self.n_vertices + self.n + self.n_omega, 
                                                            self.n_vertices + self.n + self.n_omega + n_epsilon))]
        else: epsilon_vertices = [[]]

        #print("inner_vertices: ", str(inner_vertices))
        #print("numbered_vertices: ", str(numbered_vertices))
        #print("omega_vertices: ", str(omega_vertices))
        #print("epsilon_vertices: ", str(epsilon_vertices))
        
        partition = inner_vertices + numbered_vertices + omega_vertices + epsilon_vertices
        
        #print(partition)

        return partition


    def plot_graph(self, G):
        #print("plotting final graph")
        GG = Graph(G)  # , loops=True)
        #print(self.get_partition())

        # ensure repeatable coloring ----

        partition = self.get_partition(G)
        assert len(partition) == self.n + 3

        color_dict = {}

        # inner vertices
        color_dict.update({"gray":partition[0]}) 

        # numbered hairs
        cmap = plt.get_cmap('viridis')
        for i in range(1, self.n+1): 
            assert len(partition[i]) == 1
            interpolation_value = (i - 1) / self.n
            color_dict.update({cmap(interpolation_value):partition[i]})

        # omega vertices
        color_dict.update({"red":partition[self.n + 1]}) 

        # epsilon vertices
        color_dict.update({"magenta":partition[self.n + 2]}) 
        
        #print(partition)
        #print(color_dict)

        return GG.plot(vertex_colors=color_dict, vertex_labels=True)



    def is_valid(self):
        # possible values
        l = self.genus >= 1 \
        and self.n >= 0 \
        and self.n_omega >= 0 \
        and self.n_vertices >= 0
        
        return l



    def get_generating_graphs(self):
        #print("Aggregated: get_generating_graphs")
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            print("self is not valid")
            return []
        
        # due to symmetry: 
        #   at most one eps-eps double-leg
        #   at most one eps attached to each vertex, or other label: omega or numbered
        max_epsilon = 2 + self.n_vertices + self.n + self.n_omega

        # each connected component has at least one epsilon or omega attached to it
        max_components = self.n_omega + max_epsilon

        for n_epsilon in range(max_epsilon + 1):
            for n_components in range(max_components + 1):
                for n_double_legs in range(n_components + 1):

                    V = WOHairyAggregatedGraphComplex.WOHairyAggregatedGVS(n_components=n_components,
                                         n_vertices=self.n_vertices,
                                         genus=self.genus,
                                         n=self.n,
                                         n_omega=self.n_omega,
                                         n_epsilon=n_epsilon,
                                         n_double_legs=n_double_legs,
                                         do_print=False) 

                    assert V.cohom_degree() == self.degree

                    if V.is_valid():
                        
                        """
                        print("valid-configuration: -----")
                        print("n_components", n_components)
                        print("n_vertices", self.n_vertices)
                        print("genus", self.genus)
                        print("n", self.n)
                        print("n_omega", self.n_omega)
                        print("n_epsilon", n_epsilon)
                        print("n_double_legs", n_double_legs)
                        """
                        
                        V.build_basis(ignore_existing_files=False)

                        for G in V.get_basis():
                            yield G


    def build_basis(self, progress_bar=False, ignore_existing_files=False, **kwargs):
        """Build the basis of the vector space.

        Create the basis file if the vector space is valid, otherwise skip building a basis. If there exists already
        a basis file rebuild the basis if ignore_existing_file is True, otherwise skip building a basis.
        The basis file contains a list of graph6 strings for canonically labeled graphs building a basis of the vector
        space. The canonical labeling respects the partition of the vertices.

        :param progress_bar: Option to show a progress bar (Default: False).
        :type progress_bar: bool
        :param ignore_existing_files: Option to ignore existing basis file. Ignore existing file and
                rebuild the basis if True, otherwise skip rebuilding the basis file if there exists a basis file already
                (Default: False).
        :type ignore_existing_files: bool
        :param kwargs: Accepting further keyword arguments, which have no influence.
        """
        # print("build basis ", str(self))
        if not self.is_valid():
            # Skip building a basis file if the vector space is not valid.
            return
        if (not ignore_existing_files) and self.exists_basis_file():
            # Skip building a basis file if there exists already one and ignore_existing_file is False.
            return

        generating_list = self.get_generating_graphs()

        desc = 'Build basis: ' + str(self.get_ordered_param_dict())
        # if not progress_bar:
        print(desc)
        basis_set = set()
        # for G in tqdm(generating_list, desc=desc, disable=(not progress_bar)):
        for G in generating_list:
            # For each graph G in the generating list, add the canonical labeled graph6 representation to the basis set
            # if the graph G doesn't have odd automormphisms.
            if self.get_partition(G) is None:
                autom_list = G.automorphism_group().gens()
                canonG = G.canonical_label(
                    algorithm=Parameters.canonical_label_algorithm)
            else:
                # The canonical labelling respects the partition of the vertices.
                autom_list = G.automorphism_group(
                    partition=self.get_partition(G)).gens()
                canonG = G.canonical_label(partition=self.get_partition(G), algorithm=Parameters.canonical_label_algorithm)

            canon6 = canonG.graph6_string()

            if canon6 not in basis_set:
                if not self._has_odd_automorphisms(G, autom_list):
                    basis_set.add(canon6)

        L = list(basis_set)
        L.sort()
        self._store_basis_g6(L)

    

