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
        return 'contract edges'


    # TODO: check implementation of sign!
    def operate_on(self, G):
        # self.domain is a WOHairyFinalGVS instance
        # replaces one epsilon vertex with an omega vertex
        # sem-paper: The piece δω simply removes one distinguished half-edge from the orientation

        n_epsilon = len(G.vertices()) - self.domain.n_vertices - self.domain.n - self.domain.n_omega
        assert n_epsilon >= 0

        G1 = copy(G)
        sgn = (-1)**(G.size() - n_epsilon + 1)
        image = []

        # label all edges to determine sign later
        Shared.enumerate_edges(G1)

        for i in range(n_epsilon):
            
            omega_eps_cutoff = self.domain.n_vertices + self.domain.n + self.domain.n_omega

            eps_index = omega_eps_cutoff + i
            
            # permute chosen epsilon vertex to the position where it becomes 
            # an omega vertex due to the new partition in self.domain, since:
            # self.target.n_omega = self.domain.n_omega + 1 
            
            G2 = copy(G1)

            relabeling_perm = list(range(omega_eps_cutoff)) + [eps_index] + list(range(omega_eps_cutoff, eps_index)) + list(range(eps_index+1, G1.order()))

            assert set(relabeling_perm) == set(G2.vertices())

            G2.relabel(relabeling_perm)

            sgn2 = sgn * Shared.shifted_edge_perm_sign2(G2)

            image.append((G2, sgn2))

        return image


    """
    def operate_on(self, G):
        # self is a WOHairyFinalGVS instance
        G1 = copy(G)
        sgn = (-1)**G.size()
        image = []

        # label all edges to determine sign later
        Shared.enumerate_edges(G1)

        # add one new omega vertex in position n_vertices +1
        G1.relabel(list(range(self.domain.n_vertices+1)) + list(range(self.domain.n_vertices +
                   2, self.domain.n_vertices+self.domain.n_ws+self.domain.n_hairs+2)))
        G1.add_vertex(self.domain.n_vertices+1)

        # reconnect one eps edge to the new vertex
        eps = self.domain.n_vertices
        new_w = self.domain.n_vertices+1
        # in my case we go through all epsilon vertices instead. (may be empty)
        for v in G1.neighbors(eps):
            sgn2 = sgn
            G2 = copy(G1)
            old_label = G2.edge_label(v, eps)
            G2.delete_edge(v, eps)
            G2.add_edge(v, new_w, old_label)
            sgn2 *= Shared.shifted_edge_perm_sign2(G2)
            image.append((G2, sgn2))

        return image
    """

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
    """Contract edges differential."""

    def __init__(self, sum_vector_space):
        """Initialize the contract edges differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: HairyGraphSumVS
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

    # TODO:
    # - Assert excess-compatability / degree-compatability etc after the operation
    def operate_on(self, G):

        n_old_eps = self.domain.get_n_epsilon_from_graph(G)

        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        image = []
        for (i, e) in enumerate(G.edges(labels=False,sort=True)):
            (u, v) = e

            # ensure u<v (this should be always true anyway actually)
            assert u < v

            # only edges connected to at least one internal vertex, and not connected to a numbered hair-vertex can be contracted

            # both u and v are not internal vertices
            if u >= self.domain.n_vertices: continue

            # u or v are numbered vertices        
            if u >= self.domain.n_vertices and u < self.domain.n_vertices + self.domain.n: continue
            if v >= self.domain.n_vertices and v < self.domain.n_vertices + self.domain.n: continue
        

            sgn = 1 if i % 2 == 0 else -1
            # print("sgn0",sgn)
            previous_size = G.size()

            # print("sgn1",sgn)
            G1 = copy(G)
            # label all edges to determine sign later
            Shared.enumerate_edges(G1)


            # contracting the edge (u,v)
            # we always delete the lower index vertex. This ensures that the extra vertices are never deleted

            # if both u and v are internal vertices
            if v < self.domain.n_vertices:

                # contract edge by merging vertices
                G1.merge_vertices([v, u])

                # if we have too few edges some double edges have been created => zero
                if (previous_size - G1.size()) > 1: continue

                # relabel vertices
                G1.relabel(range(len(G1.vertices())), inplace=True)
                
                # find edge permutation sign
                sgn *= Shared.shifted_edge_perm_sign2(G1)
                # print("sgn3_",sgn)
                image.append((G1, sgn))
                # image.append((Graph(G1.graph6_string()), sgn))
                # print("hmm0:", G.graph6_string(), G1.graph6_string())


            # if u is an internal vertex and v is a omega-vertex
            # the second vertex is now an omega-vertex, so we need to merge the vertex with the eps vertex
            # after reconnecting one of the edges to omega
            # we assume that u != eps, because eps-omega-edges cannot be contracted
            elif u < self.domain.n_vertices \
                and v >= self.domain.n_vertices + self.domain.n \
                and v < self.domain.n_vertices + self.domain.n + self.domain.n_omega:

                G1.delete_edge(u, v)

                # loop over neighbors w of u to be connected to omega
                # TODO: Question: if we delete all edges adjacent to u, 
                #       will it not just be an isolated vertex?
                for w in G1.neighbors(u):
                    G2 = copy(G1)
                    sgn2 = sgn

                    # reconnect the w-v-edge to omega (i.e., to v)
                    old_label = G2.edge_label(u, w)
                    G2.delete_edge(u, w)
                    G2.add_edge(w, v, old_label)

                    # we now need to split u into separate epsilon-vertices
                    # for each remaining adjecent vertex s of u
                    n_new_eps = len(G2.neighbors(u))
                    for (j, s) in enumerate(G2.neighbors(u)): 
                        G2.delete_edge(u, s)

                        assert set(G2.vertices()) == set(range(len(G2.vertices())))
                        # since we add the new vertex at the end, it automatcally corresponds to epsilon
                        new_vertex_label = len(G2.vertices()) 

                        old_size = G2.order()
                        #print(G2.vertices())
                        G2.add_vertex(new_vertex_label)
                        assert G2.order() == old_size + 1, (G2.order(), old_size)

                        label = 0 #TODO!
                        G2.add_edge(new_vertex_label, s, label)

                    G2.delete_vertex(u)
                        

                    # in case we have too few edges some double edges have been created => zero
                    if (previous_size - G2.size()) > 1: continue

                    G2.relabel(range(len(G2.vertices())), inplace=True)
                    # find edge permutation sign
                    sgn2 *= Shared.shifted_edge_perm_sign2(G2)

                    # sanity checks
                    n_eps_in_target = self.target.get_n_epsilon_from_graph(G2)
                    assert n_eps_in_target == n_old_eps + n_new_eps
                    assert G2.order() == self.target.n_vertices + self.target.n + self.target.n_omega + n_eps_in_target


                    image.append((G2, sgn2))


        return image
    
    """
    def operate_on(self, G):
        # print("operate on:", G.graph6_string(),
            #   self.domain.get_ordered_param_dict())
        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        image = []
        for (i, e) in enumerate(G.edges(labels=False,sort=True)):
            (u, v) = e

            # ensure u<v (this should be always true anyway actually)
            if u > v:
                u, v = v, u

            # only edges connected to at least one internal vertex, and not connected to a numbered hair-vertex can be contracted
            if u >= self.domain.n_vertices or v >= self.domain.n_vertices+self.domain.n_ws+1:
                continue

            sgn = 1 if i % 2 == 0 else -1
            # print("sgn0",sgn)
            previous_size = G.size()
            previous_has_tadpole = (
                previous_size - self.domain.n_vertices - self.domain.n_hairs < self.domain.n_loops)
            sgn *= -1 if previous_has_tadpole else 1
            # print("sgn1",sgn)
            G1 = copy(G)
            # label all edges to determine sign later
            Shared.enumerate_edges(G1)

            # we always delete the lower index vertex. This ensures that the extra vertices are never deleted
            if v <= self.domain.n_vertices:
                G1.merge_vertices([v, u])
                if (previous_size - G1.size()) != 1:
                    continue
                G1.relabel(range(0, self.domain.n_vertices+self.domain.n_ws +
                           self.domain.n_hairs), inplace=True)
                # find edge permutation sign
                sgn *= Shared.shifted_edge_perm_sign2(G1)
                # print("sgn3_",sgn)
                image.append((G1, sgn))
                # image.append((Graph(G1.graph6_string()), sgn))
                # print("hmm0:", G.graph6_string(), G1.graph6_string())
            elif u < self.domain.n_vertices and v >= self.domain.n_vertices+1:
                # the second vertex is now an omega-vertex, so we need to merge the vertex with the eps vertex
                # after reconnecting one of the edges to omega
                # we assume that u != eps, because eps-omega-edges cannot be contracted
                G1.delete_edge(u, v)
                # special care must be taken since a tadpole could be created at eps
                # and this is true iff there is an edge u-eps
                eps = self.domain.n_vertices
                # new_has_tadpole = G1.has_edge(u, eps)
                # # double tadpole => zero
                # if new_has_tadpole and previous_has_tadpole:
                #     continue
                # if new_has_tadpole:
                #     # remove the edge and compute the appropriate sign
                #     k = G1.edge_label(u, eps)
                #     G1.delete_edge(u, eps)
                #     sgn *= 1 if ((k % 2 == 0) == (k < i)) else -1

                # loop over neighbors w to be connected to omega
                for w in G1.neighbors(u):
                    G2 = copy(G1)
                    sgn2 = sgn
                    # reconnect the w-v-edge to omega (i.e., to v)
                    old_label = G2.edge_label(u, w)
                    G2.delete_edge(u, w)
                    G2.add_edge(w, v, old_label)

                    # we want to merge u and eps... however, this might create a tadpole
                    new_has_tadpole = G2.has_edge(u, eps)
                    if new_has_tadpole and previous_has_tadpole:
                        continue
                    if new_has_tadpole:
                        # remove the edge and compute the appropriate sign
                        k = G2.edge_label(u, eps)
                        G2.delete_edge(u, eps)
                        sgn2 *= 1 if ((k % 2 == 0) == (k < i)) else -1

                    # now merge u and eps
                    G2.merge_vertices([eps, u])
                    # in case we have too few edges some double edges have been created => zero
                    if (previous_size - G2.size()) != (2 if new_has_tadpole else 1):
                        continue
                    G2.relabel(range(0, self.domain.n_vertices +
                               self.domain.n_hairs+self.domain.n_ws), inplace=True)
                    # find edge permutation sign
                    sgn2 *= Shared.shifted_edge_perm_sign2(G2)
                    # sanity checks
                    if G2.order() != self.target.n_vertices+self.target.n_hairs+self.target.n_ws+1:
                        print("Error contract:", G.graph6_string(),
                              G2.graph6_string())
                    # else:
                    #     print("hmm:", G.graph6_string(), G2.graph6_string())
                    image.append((G2, sgn2))

        return image
    """
    
    
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
        """Initialize the graph complex.

        """
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

    def print_dim_and_eulerchar(self):
        for n_omega in self.omega_range:
            for n in self.n_range:
                for genus in self.genus_range:
                    ds = [WOHairyFinalGVS(genus=genus, n=n, n_omega=n_omega, degree=degree).get_dimension()
                          for degree in self.degree_range]
                    eul = sum([(1 if j % 2 == 0 else -1) *
                              d for j, d in enumerate(ds)])
                    print("Dimensions (n_omega,n,genus) ", n_omega,
                          n, genus, ":", ds, "Euler", eul)

    # TODO: Insert contribution from EpsToOmegaGo
    def print_cohomology_dim(self):
        for n_omega in self.omega_range:
            for n in self.n_range:
                for genus in self.genus_range:
                    cohomdict = {}
                    for degree in self.degree_range:
                        D1 = ContractEdgesGO.generate_operator(degree, genus, n, n_omega)
                        D2 = ContractEdgesGO.generate_operator(degree+1, genus, n, n_omega)
                        try:
                            d = WOHairyFinalGVS(genus=genus, n=n, n_omega=n_omega, degree=degree).get_dimension()
                            r1 = D1.get_matrix_rank()
                            r2 = D2.get_matrix_rank()
                            cohomdict[degree] = d-r1-r2
                        except:
                            pass

                    print("Cohomology Dimensions (w,n,genus) ",
                          n_omega, n, genus, ":", cohomdict)


