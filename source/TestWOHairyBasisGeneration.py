
# testing of basis-generation functionality using simple elementary examples and euler-characteristic
# for testing to be fully rigurous, the "data"-folder should be deleted before running the tests


import unittest
from sage.all import StandardTableaux
from WOHairyGC_Pascal import WOHairyComponentGVS, WOHairyAggregatedGVS, WOHairyGVS


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

        V = WOHairyAggregatedGVS(n_components=n_components, n_vertices=n_vertices, genus=genus, n=n, n_omega=n_omega, n_epsilon=n_epsilon, n_double_legs=n_double_legs)

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
                                                    n_double_legs=n_double_legs)
                        
                        assert V.excess == excess

                        V.build_basis(ignore_existing_files=True)

                        assert V.get_dimension() == 1, (excess, n_omega, n)





class TestFinalGeneration(unittest.TestCase):

    @staticmethod
    def test_basis_len(test_name, test_basis_len, genus, n, n_omega, degree):

        V = WOHairyGVS(genus=genus, n=n, n_omega=n_omega, degree=degree)

        if test_basis_len > 0: assert V.is_valid(), test_name

        V.build_basis(ignore_existing_files=True)

        assert V.get_dimension() == test_basis_len, test_name


    def test_bases_by_excess(self):

        # testing parameters (g, n, degree) in excess 0 with notrivial bases
        self.test_basis_len("excess 0: B_1_11", 1, 1, 11, 11, 11)
        self.test_basis_len("excess 0: B_3_8", 1, 3, 8, 11, 14)
        self.test_basis_len("excess 0: B_5_5", 1, 5, 5, 11, 17)
        self.test_basis_len("excess 0: B_7_2", 1, 7, 2, 11, 20)

        # excess 1
        self.test_basis_len("excess 1: B_2_10", 10, 2, 10, 11, 13)
        self.test_basis_len("excess 1: B_6_4", 5, 6, 4, 11, 19)

        # excess 2
        self.test_basis_len("excess 2: B_7_3", 16, 7, 3, 11, 20)





class TestEulerCharcteristic(unittest.TestCase):

    
    

    @staticmethod
    def test_euler_char(genus, n, coefficients, diagrams, excess):

        print("testing for (g,n)=", (genus, n), "-----------------------")

        computed_excess = 3*genus + 2*n - 25

        if computed_excess < 0: assert excess < 0
        else: assert computed_excess == excess

        euler_char = WOHairyGVS.compute_euler_char(genus, n)
        
        assert len(coefficients) == len(diagrams)

        euler_char_reference = 0
        for coefficient, diagram in zip(coefficients, diagrams):
            euler_char_reference += coefficient * StandardTableaux(diagram).cardinality()

        print("euler_char:", euler_char)
        print("euler_char_reference:", euler_char_reference)
        assert euler_char == euler_char_reference, (euler_char, euler_char_reference)


    def test_eulerChar(self):
        # excess < 0
        for i in range(5): self.test_euler_char(5, i, coefficients=[], diagrams=[], excess=-1)
        for i in range(4): self.test_euler_char(6, i, coefficients=[], diagrams=[], excess=-1)
        for i in range(2): self.test_euler_char(7, i, coefficients=[], diagrams=[], excess=-1)
        self.test_euler_char(8, 0, coefficients=[], diagrams=[], excess=-1)

        # excess = 0
        self.test_euler_char(5, 5, coefficients=[-1], diagrams=[[1,1,1,1,1]], excess=0)
        self.test_euler_char(7, 2, coefficients=[1], diagrams=[[1,1]], excess=0)

        # excess = 1
        self.test_euler_char(6, 4, coefficients=[-1], diagrams=[[2,1,1]], excess=1)
        self.test_euler_char(8, 1, coefficients=[], diagrams=[], excess=1)

        # excess = 2
        self.test_euler_char(5, 6, coefficients=[-1, 2], diagrams=[[1,1,1,1,1,1], [3,1,1,1]], excess=2)
        self.test_euler_char(7, 3, coefficients=[1, -2], diagrams=[[1,1,1], [3]], excess=2)
        self.test_euler_char(9, 0, coefficients=[1], diagrams=[[]], excess=2)

        # excess = 3
        self.test_euler_char(8, 2, coefficients=[1], diagrams=[[2]], excess=3)
        self.test_euler_char(6, 5, coefficients=[-1, 1, 1, 2], diagrams=[[2,1,1,1], [3,1,1], [3,2], [4,1]], excess=3)

        # excess = 4
        self.test_euler_char(9, 1, coefficients=[1], diagrams=[[1]], excess=4)
        self.test_euler_char(7, 4, coefficients=[2, -2], diagrams=[[2,2], [3,1]], excess=4)


        # excess = 5
        self.test_euler_char(10, 0, coefficients=[-2], diagrams=[[]], excess=5)
        self.test_euler_char(8, 3, coefficients=[4], diagrams=[[1,1,1]], excess=5)
        #self.test_euler_char(6, 6, coefficients=[-4, -2, -3, -3, 1, -2, -2, 2, -3, -2 ,-3], 
        #                diagrams=[[1,1,1,1,1,1], [2,1,1,1,1], [2,2,1,1], [2,2,2], 
        #                        [3,1,1,1], [3,2,1], [3,3], [4,1,1], [4,2], [5,1], [6]], excess=5)

        # excess = 6
        self.test_euler_char(9, 2, coefficients=[4], diagrams=[[2]], excess=6)
        #self.test_euler_char(7, 5, coefficients=[-1, -7, -6, 2, 1, 4], diagrams=[[1,1,1,1,1], [2,1,1,1], [3,1,1], [3,2], [4,1], [5]], excess=6)

        # excess = 7
        self.test_euler_char(10, 1, coefficients=[-6], diagrams=[[1]], excess=7)
        #self.test_euler_char(8, 4, coefficients=[10, 2, -2, -10, -3], diagrams=[[1,1,1,1], [2,1,1], [2,2], [3,1], [4]], excess=7)

        # excess = 8
        self.test_euler_char(11, 0, coefficients=[2], diagrams=[[]], excess=8)
        #self.test_euler_char(9, 3, coefficients=[1, 11, 4], diagrams=[[1,1,1], [2,1], [3]], excess=8)

        # excess = 9
        #self.test_euler_char(10, 2, coefficients=[-9, -2], diagrams=[[1,1], [2]], excess=9)

        # excess = 10
        #self.test_euler_char(11, 1, coefficients=[1], diagrams=[[1]], excess=10)
        




def suite():
    suite = unittest.TestSuite()

    suite.addTest(TestComponentGeneration('test_double_leg_variants_WOHairyComponentGVS'))
    suite.addTest(TestComponentGeneration('test_single_vertex_variants_WOHairyComponentGVS'))
    suite.addTest(TestComponentGeneration('test_trees_WOHairyComponentGVS'))

    suite.addTest(TestComponentAggregation('test_double_leg_variants_WOHairyAggregatedGVS'))
    suite.addTest(TestComponentAggregation('test_single_vertex_variants_WOHairyAggregatedGVS'))
    suite.addTest(TestComponentAggregation('test_single_mutliple_components_WOHairyAggregatedGVS'))
    suite.addTest(TestComponentAggregation('test_trees_WOHairyAggregatedGVS'))

    suite.addTest(TestFinalGeneration('test_bases_by_excess'))

    suite.addTest(TestEulerCharcteristic('test_eulerChar'))

    return suite


if __name__ == '__main__':
    print("\n#######################################\n" + "----- Start test suite WOHairy Baisis Generation -----")
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
