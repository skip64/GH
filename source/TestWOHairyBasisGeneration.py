import unittest

# testing of graph generation functionality using simple elementary examples 
# and euler-characteristic

from WOHairyComponentGraphComplex import WOHairyComponentGVS
from WOHairyAggregatedGraphComplex import WOHairyAggregatedGVS
from WOHairyFinalGraphComplex import WOHairyFinalGVS





class TestComponentGeneration(unittest.TestCase):


    @staticmethod
    def test_basis_len(test_name, basis_len, n_vertices, n_loops, n, n_omega, n_epsilon):

        V = WOHairyComponentGVS(n_vertices=n_vertices, n_loops=n_loops, n=n, n_omega=n_omega, n_epsilon=n_epsilon)

        if basis_len > 0: assert V.is_valid(), test_name

        V.build_basis(ignore_existing_files=True)

        assert V.get_dimension() == basis_len, test_name


    def test_double_leg_variants_WOHairyComponentGVS(self):
        # testing double-leg configurations
        self.test_basis_len("omega-omega", basis_len=0, n_vertices=0, n_loops=0, n=0, n_omega=2, n_epsilon=0)
        self.test_basis_len("epsilon-epsilon", basis_len=1, n_vertices=0, n_loops=0, n=0, n_omega=0, n_epsilon=2)
        self.test_basis_len("omega-epsilon", basis_len=1, n_vertices=0, n_loops=0, n=0, n_omega=1, n_epsilon=1)
        self.test_basis_len("leg-omega", basis_len=1, n_vertices=0, n_loops=0, n=1, n_omega=1, n_epsilon=0)
        self.test_basis_len("leg-epsilon", basis_len=1, n_vertices=0, n_loops=0, n=1, n_omega=0, n_epsilon=1)


    def test_single_vertex_variants_WOHairyComponentGVS(self):
        # testing graphs with one vertex
        self.test_basis_len("eps-omega-eps", basis_len=0, n_vertices=1, n_loops=0, n=0, n_omega=1, n_epsilon=2)
        self.test_basis_len("3-omega", basis_len=1, n_vertices=1, n_loops=0, n=0, n_omega=3, n_epsilon=0)
        self.test_basis_len("3-epsilon", basis_len=0, n_vertices=1, n_loops=0, n=0, n_omega=0, n_epsilon=3)


    def test_trees_WOHairyComponentGVS(self):
        # testing trees built from omegas and numbered legs on a single vertex

        n_components = 1
        n_vertices = 1
        n_epsilon = 0
        n_double_legs = 0

        for excess in range(6):

            for n_omega in range(excess + 4):
                
                genus = n_omega

                # excess = n_omega - 3 + 2*n -> 2n = excess - n_omega + 3
                two_n = excess - n_omega + 3

                if two_n % 2 == 0: 
                    n = int(two_n / 2)

                    # computing loop-order
                    n_hairs = n + n_omega + n_epsilon
                    n_edges = genus + n_vertices + n + n_double_legs - n_hairs - 1
                    n_loops = n_edges - n_vertices + (n_components - n_double_legs)

                    if n_omega + n >= 3 and n_omega >= 1:

                        V = WOHairyComponentGVS(n_vertices=n_vertices, 
                                                n_loops=n_loops, 
                                                n=n, n_omega=n_omega, n_epsilon=n_epsilon)

                        V.build_basis(ignore_existing_files=True)

                        assert V.get_dimension() == 1, (excess, n_omega, n)



class TestComponentAggregation(unittest.TestCase):


    @staticmethod
    def test_basis_len(test_name, basis_len, inner_edges, n_components, n_vertices, genus, n, n_omega, n_epsilon, n_double_legs):

        V = WOHairyAggregatedGVS(n_components=n_components, n_vertices=n_vertices, genus=genus, n=n, n_omega=n_omega, n_epsilon=n_epsilon, n_double_legs=n_double_legs, do_print=False)

        assert V.n_edges == inner_edges, inner_edges

        if basis_len > 0: assert V.is_valid(), test_name

        V.build_basis(ignore_existing_files=True)

        assert V.get_dimension() == basis_len, test_name


    def test_double_leg_variants_WOHairyAggregatedGVS(self):
        # testing double-leg configurations

        n_components = 1
        n_vertices = 0
        n_double_legs = 1

        self.test_basis_len("omega-omega", basis_len=0, inner_edges=0, genus=2, n=0, n_omega=2, n_epsilon=0, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)
        self.test_basis_len("epsilon-epsilon", basis_len=1, inner_edges=0, genus=2, n=0, n_omega=0, n_epsilon=2, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)
        self.test_basis_len("omega-epsilon", basis_len=1, inner_edges=0, genus=2, n=0, n_omega=1, n_epsilon=1, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)
        self.test_basis_len("leg-omega", basis_len=1, inner_edges=0, genus=1, n=1, n_omega=1, n_epsilon=0, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)
        self.test_basis_len("leg-epsilon", basis_len=1, inner_edges=0, genus=1, n=1, n_omega=0, n_epsilon=1, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)


    def test_single_vertex_variants_WOHairyAggregatedGVS(self):
        # testing graphs with one vertex

        n_components = 1
        n_vertices = 1
        n_double_legs = 0

        self.test_basis_len("eps-omega-eps", basis_len=0, inner_edges=0, genus=3, n=0, n_omega=1, n_epsilon=2, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)
        self.test_basis_len("3-omega", basis_len=1, inner_edges=0, genus=3, n=0, n_omega=3, n_epsilon=0, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)
        self.test_basis_len("3-epsilon", basis_len=0, inner_edges=0, genus=3, n=0, n_omega=0, n_epsilon=3, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)


    def test_single_mutliple_components_WOHairyAggregatedGVS(self):
        # testing graphs with multiple components ---

        self.test_basis_len("omega-1 & omega-2", basis_len=1, inner_edges=0, n_components=2, n_vertices=0, genus=1, n=2, n_omega=2, n_epsilon=0, n_double_legs=2)
        self.test_basis_len("omega-1 & omega-2 & 2x3-omega", basis_len=1, inner_edges=0, n_components=4, n_vertices=2, genus=5, n=2, n_omega=8, n_epsilon=0, n_double_legs=2)



    def test_trees_WOHairyAggregatedGVS(self):
        # testing trees built from omegas and numbered legs on a single vertex

        n_components = 1
        n_vertices = 1
        n_epsilon = 0
        n_double_legs = 0

        for excess in range(10):

            for n_omega in range(excess + 4):

                genus = n_omega

                # excess = n_omega - 3 + 2*n -> 2n = excess - n_omega + 3
                two_n = excess - n_omega + 3

                if two_n % 2 == 0: 
                    n = int(two_n / 2)

                    if n_omega + n >= 3 and n_omega >= 1:

                        V = WOHairyAggregatedGVS(n_components=n_components, 
                                                    n_vertices=n_vertices, 
                                                    genus=genus, 
                                                    n=n, n_omega=n_omega, n_epsilon=n_epsilon,
                                                    n_double_legs=n_double_legs,
                                                    do_print=False)
                        
                        assert V.excess == excess

                        V.build_basis(ignore_existing_files=True)

                        assert V.get_dimension() == 1, (excess, n_omega, n)


class TestFinalGeneration(unittest.TestCase):

    def setUp(self):
        pass


class TestEulerCharcteristic(unittest.TestCase):

    def setUp(self):
        pass


def suite():
    suite = unittest.TestSuite()

    suite.addTest(TestComponentGeneration('test_double_leg_variants_WOHairyComponentGVS'))
    suite.addTest(TestComponentGeneration('test_single_vertex_variants_WOHairyComponentGVS'))
    suite.addTest(TestComponentGeneration('test_trees_WOHairyComponentGVS'))

    suite.addTest(TestComponentAggregation('test_double_leg_variants_WOHairyAggregatedGVS'))
    suite.addTest(TestComponentAggregation('test_single_vertex_variants_WOHairyAggregatedGVS'))
    suite.addTest(TestComponentAggregation('test_single_mutliple_components_WOHairyAggregatedGVS'))
    suite.addTest(TestComponentAggregation('test_trees_WOHairyAggregatedGVS'))


    return suite


if __name__ == '__main__':
    print("\n#######################################\n" + "----- Start test suite WOHairy Baisis Generation -----")
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
