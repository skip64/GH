import SymmetricGraphComplex
import GraphOperator
import os
import Parameters
import StoreLoad
import math
import Shared
from copy import copy
from WOHairyBasisGeneration import WOHairyFinalGVS
import GraphVectorSpace
import itertools
import GraphComplex

graph_type = "wohairy_final"


"""
Implementation-Notes:
- do not need sub_type: (even / odd) -> will try to remove...
"""



class WOHairyGraphSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of graph vector spaces with specified number of omega hairs.
    """

    def __init__(self, genus_range, n_range, n_omega_range, degree_range):
        """Initialize the sum vector space.
        """
        self.genus_range = genus_range
        self.n_range = n_range
        self.n_omega_range = n_omega_range
        self.degree_range = degree_range

        vs_list = [WOHairyFinalGVS(genus, n, n_omega, degree) for
                   (genus, n, n_omega, degree) in itertools.product(self.genus_range, self.n_range, self.n_omega_range, self.degree_range)]
        super(WOHairyGraphSumVS, self).__init__(vs_list)

    def get_type(self):
        return 'wohairy graphs'

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('genus', self.genus_range), ('n', self.n_range), ('n_omega', self.n_omega_range), ('degree', self.degree_range)])

    def get_info_plot_path(self):
        s = "info_vector_space_%s" % graph_type
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)





# ------- Operators --------

class EpsToOmegaGO(SymmetricGraphComplex.SymmetricGraphOperator):
    """Operator that makes one eps into an omega hair.
    """

    def __init__(self, domain, target):
        super(EpsToOmegaGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        return (domain.genus == target.genus 
                and domain.n == target.n
                and domain.n_omega + 1 == target.n_omega 
                and domain.degree - 1 == target.degree)

    @classmethod
    def generate_operator(cls, genus, n, n_omega, degree):
        domain = WOHairyFinalGVS(genus, n, n_omega, degree)
        target = WOHairyFinalGVS(genus, n, n_omega + 1, degree - 1)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "epstowD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_rank_file_path(self):
        s = "epstowD%d_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_type(self):
        return 'eps to omega'


    # TODO: check implementation of sign!
    def operate_on(self, G):
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
        notation: "i:0":
        - "i" denotes the "physical" vertex in the graph
        - "0" denotes the corresponding vertex-label
        1.  G1.vertices():  i:0, o1:1, o2:2, e1:3, e2:4, e3:5, e4:6, e5:7
        2.  eps_index = 5 -> eps_vertex e3
        3.  we now want to swap vertex e3 with e1: 
            relabeling_perm: [0, 1, 2, 5, 4, 3, 6, 7] 
            G2.vertices():  i:0, o1:1, o2:2, e1:5, e2:4, e3:3, e4:6, e5:7
        4.  Now, since each hair-vertex (o1, o2, e1, e2, e3, e4, e5) corresponds to exactly one edge (the hair),
            the induced edge permutation is simply a swap of (i-e3) with (i-e1)
            -> sign = -1
        """

        n_epsilon = len(G.vertices()) - self.domain.n_vertices - self.domain.n - self.domain.n_omega
        assert n_epsilon >= 0

        G1 = copy(G)
        sgn = 1
        image = []

        # label all edges to determine sign later
        Shared.enumerate_edges(G1)

        for i in range(n_epsilon):
            
            omega_eps_cutoff = self.domain.n_vertices + self.domain.n + self.domain.n_omega
            eps_index = omega_eps_cutoff + i
            
            # permute chosen epsilon vertex to the position where it becomes an omega vertex due to the new partition in self.domain
            # if i == 0: then omega_eps_cutoff == eps_index
            G2 = copy(G1)
            if i > 0:
                # relabeling_perm swaps omega_eps_cutoff and eps_index
                relabeling_perm = list(range(omega_eps_cutoff)) + [eps_index] + list(range(omega_eps_cutoff+1, eps_index)) + [omega_eps_cutoff] + list(range(eps_index+1, G1.order()))
                #print(relabeling_perm)
                assert set(relabeling_perm) == set(G2.vertices())
                G2.relabel(relabeling_perm)

            # sign equals the perm_sign of the edge-permutation induced by the relabeling of the vertices
            sgn = Shared.shifted_edge_perm_sign2(G2)
            #print("sgn: ", sgn)
            image.append((G2, sgn))

        return image


    def get_work_estimate(self):
        return 0

    def restrict_to_isotypical_component(self, rep_index):
        return RestrictedEpsToOmegaGO(self, rep_index)



class RestrictedEpsToOmegaGO(SymmetricGraphComplex.SymmetricRestrictedOperatorMatrix):
    # def __init__(opD, opP):

    def get_matrix_file_path(self):
        s = "contractD%d_%d_%d_r%d.txt" % (
            self.domain.vs.get_ordered_param_dict().get_value_tuple() + (self.rep_index,))
        return os.path.join(Parameters.data_dir, graph_type, self.opD.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD%d_%d_%d_r%d_rank.txt" % (
            self.domain.vs.get_ordered_param_dict().get_value_tuple() + (self.rep_index,))
        return os.path.join(Parameters.data_dir, graph_type, self.opD.sub_type, s)

    def get_work_estimate(self):
        return self.opD.get_work_estimate()

    def is_match(self, domain, target):
        return EpsToOmegaGO.is_match(domain.vs, target.vs) and domain.rep_index == target.rep_index


class EpsToOmegaD(GraphOperator.Differential):

    def __init__(self, sum_vector_space):
        """Initialize the eps to omega differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: 
        """
        super(EpsToOmegaD, self).__init__(sum_vector_space,
                                          EpsToOmegaGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'EpsToOmega'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_epstoomega_D_%s" % (graph_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_cohomology_web_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_epstoomega_D_%s" % (graph_type)
        return os.path.join(Parameters.web_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_epstoomega_D_%s" % (graph_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)
    



class RestrictedEpsToOmegaD(SymmetricGraphComplex.SymmetricDifferential):

    def get_type(self):
        return 'isotypical epstoomega'

    def get_cohomology_plot_path(self):
        sub_type = self.diff.sum_vector_space.sub_type
        s = "cohomology_dim_epstoomega_D_iso_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.diff.sum_vector_space.sub_type
        s = "info_epstoomega_D_iso_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)
    




class ContractEdgesGO(SymmetricGraphComplex.SymmetricGraphOperator):
    """Contract edges graph operator.
    """

    def __init__(self, domain, target):
        super(ContractEdgesGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        return (domain.genus == target.genus 
                and domain.n == target.n
                and domain.n_omega == target.n_omega 
                and domain.degree - 1 == target.degree)

    @classmethod
    def generate_operator(cls, genus, n, n_omega, degree):
        domain = WOHairyFinalGVS(genus, n, n_omega, degree)
        target = WOHairyFinalGVS(genus, n, n_omega, degree-1)
        return cls(domain, target)


    def get_matrix_file_path(self):
        s = "contractD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_rank_file_path(self):
        s = "contractD%d_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD%d_%d_%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_type(self):
        return 'contract edges'


    def operate_on(self, G):
        
        n_vertices = self.domain.n_vertices
        n = self.domain.n
        n_omega = self.domain.n_omega
        n_epsilon = self.domain.get_n_epsilon_from_graph(G)

        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        image = []
        for (i, e) in enumerate(G.edges(labels=False, sort=True)):

            #print("contracting edge (u,v) =", e, "###########")
            
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

            # label all edges to determine sign later
            G1 = copy(G)
            Shared.enumerate_edges(G1)


            # contracting the edge (u,v)
            # we always delete the lower index vertex. This ensures that the extra vertices are never deleted

            # if both u and v are internal vertices
            if v < n_vertices:
                
                #print("CASE: inner_vertex - inner_vertex")

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
                image.append((G1, sgn))
                # image.append((Graph(G1.graph6_string()), sgn))
                # print("hmm0:", G.graph6_string(), G1.graph6_string())


            # if u is an internal vertex and v is an epsilon-vertex
            elif u < n_vertices \
                and v >= n_vertices + n + n_omega:
                
                #print("CASE: inner_vertex - epsilon")

                G1.delete_vertex(v)

                # split u into separate epsilon-vertices
                for w in G1.neighbors(u):

                    #print("picking neighbour w =", w, " ---")
                    
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
                image.append((G1, sgn))


            # if u is an internal vertex and v is a omega-vertex
            # the second vertex is now an omega-vertex, so we need to merge the vertex with the eps vertex
            # after reconnecting one of the edges to omega
            # we assume that u != eps, because eps-omega-edges cannot be contracted
            elif u < n_vertices \
                and v >= n_vertices + n \
                and v < n_vertices + n + n_omega:
                
                #print("CASE: inner_vertex - omega")

                G1.delete_vertex(v)

                # pick vertex w which will be connected to omega
                for (j, w) in enumerate(G1.neighbors(u)):

                    #print("picking neighbour w =", w, " ---")
                    G2 = copy(G1)
                    sgn2 = sgn

                    u_w_label = G2.edge_label(u, w)
                    G2.delete_edge(u, w)

                    # why v is convenient:
                    # - v was the label of some omega-vertex 
                    #   ->  it will again correspond to an omega vertex in the new partition in self.target
                    #       as long as we insure to add the new epsilons after it
                    # - v has been deleted and is hence a "free" vertex-label
                    new_omega_label = v
                    assert not new_omega_label in G2.vertices()

                    G2.add_vertex(new_omega_label)
                    G2.add_edge(w, new_omega_label, u_w_label)

                    # all other vertices s will be connected to epsilons
                    n_new_eps = len(G2.neighbors(u))
                    for s in G2.neighbors(u): 
                        
                        #print("picking neighbour s =", s, " ---")
                    
                        u_s_label = G2.edge_label(u, s)
                        G2.delete_edge(u, s)

                        new_eps_label = max(G2.vertices()) + 1
                        assert not new_eps_label in G2.vertices()

                        G2.add_vertex(new_eps_label)
                        G2.add_edge(s, new_eps_label, u_s_label)
                    
                    G2.delete_vertex(u) 

                    # relabel vertices back to 0,1,2,...,k
                    G2.relabel(range(len(G2.vertices())), inplace=True)
                    
                    sgn2 *= Shared.shifted_edge_perm_sign2(G2)

                    # sanity-checks
                    n_eps_in_target = self.target.get_n_epsilon_from_graph(G2)
                    assert n_eps_in_target == n_epsilon + n_new_eps
                    assert G2.order() == self.target.n_vertices + self.target.n + self.target.n_omega + n_eps_in_target

                    image.append((G2, sgn2))


        return image
    
    
    
    def get_work_estimate(self):
        # TODO
        return 0

    def restrict_to_isotypical_component(self, rep_index):
        return RestrictedContractEdgesGO(self, rep_index)


class RestrictedContractEdgesGO(SymmetricGraphComplex.SymmetricRestrictedOperatorMatrix):
    # def __init__(opD, opP):

    def get_matrix_file_path(self):
        s = "contractD%d_%d_%d_r%d.txt" % (
            self.domain.vs.get_ordered_param_dict().get_value_tuple() + (self.rep_index,))
        return os.path.join(Parameters.data_dir, graph_type, self.opD.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD%d_%d_%d_r%d_rank.txt" % (
            self.domain.vs.get_ordered_param_dict().get_value_tuple() + (self.rep_index,))
        return os.path.join(Parameters.data_dir, graph_type, self.opD.sub_type, s)

    def get_work_estimate(self):
        return self.opD.get_work_estimate()

    def is_match(self, domain, target):
        return ContractEdgesGO.is_match(domain.vs, target.vs) and domain.rep_index == target.rep_index


class ContractEdgesD(GraphOperator.Differential):
    """Contract edges differential."""

    def __init__(self, sum_vector_space):
        """Initialize the contract edges differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: HairyGraphSumVS
        """
        super(ContractEdgesD, self).__init__(sum_vector_space,
                                             ContractEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'contract edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_D_%s" % (graph_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_cohomology_web_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_D_%s" % (graph_type)
        return os.path.join(Parameters.web_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_contract_D_%s" % (graph_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)




class RestrictedContractEdgesD(SymmetricGraphComplex.SymmetricDifferential):

    def get_type(self):
        return 'isotypical contract edges'

    def get_cohomology_plot_path(self):
        sub_type = self.diff.sum_vector_space.sub_type
        s = "cohomology_dim_contract_D_iso_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.diff.sum_vector_space.sub_type
        s = "info_contract_D_iso_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)





class SymmProjector(SymmetricGraphComplex.SymmetricProjectionOperator):
    """This class encodes the projector to an isotypical component of the symmetric group action
        by permuting numbered hairs.
        Warning: The matrix stores not the projector, but projector * n_hairs! / rep_dimension??, to have integral matrices.

    Attributes:
        - sub_type(str): Graphs sub type of the domain.
    """

    def __init__(self, domain, rep_index):
        """Initialize the domain and target vector space of the contract edges graph operator.

        : param domain: Domain vector space of the operator.
        : type domain: HairyGraphVS
        : param rep_index: The index of the representation in the list produced by Partitions(h).
        : type rep_index: int
        """
        self.sub_type = domain.sub_type

        super(SymmProjector, self).__init__(domain, rep_index)

    def get_ordered_param_dict2(self):
        do = self.domain
        return Shared.OrderedDict([('vertices', do.n_vertices), ('loops', do.n_loops), ('hairs', do.n_hairs), ('ws', do.n_ws), ('rep_index', self.rep_index)])

    def get_matrix_file_path(self):
        s = "projectionO%d_%d_%d_%d_%d.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "projectionO%d_%d_%d_%d_%d_rank.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "projectionO%d_%d_%d_%d_%d.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "projectionO%d_%d_%d_%d_%d.txt.rank.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)


# ------- Graph Complex --------
class WOHairyGC(GraphComplex.GraphComplex):
    """Graph complex for hairy graphs.
    """

    def __init__(self, genus_range, n_range, omega_range, degree_range, differentials):
        """Initialize the graph complex."""
        self.genus_range = genus_range
        self.n_range = n_range
        self.omega_range = omega_range
        self.degree_range = degree_range

        sum_vector_space = WOHairyGraphSumVS(
            self.genus_range, self.n_range, self.omega_range, self.degree_range)
        
        differential_list = []

        if not set(differentials).issubset(['contract', 'contract_iso', 'epstoomega', 'epstoomega_iso']):
            raise ValueError(
                "Differentials for hairy graph complex: 'contract'")
        
        contract_edges_dif = ContractEdgesD(sum_vector_space)
        epstoomega_dif = EpsToOmegaD(sum_vector_space)

        if 'contract' in differentials:
            differential_list.append(contract_edges_dif)

        if 'contract_iso' in differentials:
            contract_iso_edges_dif = RestrictedContractEdgesD(
                contract_edges_dif)
            differential_list.append(contract_iso_edges_dif)
            print("Attention: contract_iso operates on nonzero cohomology entries only, so they need to be computed before!")

        if 'epstoomega' in differentials:
            differential_list.append(epstoomega_dif)
        super(WOHairyGC, self).__init__(sum_vector_space, differential_list)

    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))



    



