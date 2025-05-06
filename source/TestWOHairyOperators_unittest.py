# testing of operator-functionality
# for testing to be fully rigurous, the "data"-folder should be deleted before running the tests


import unittest
from WOHairyGC_Pascal import WOHairyGC
from sage.all import *
from sage.all import StandardTableaux


# TODO: put into WOHairyGC_Pascal with assertions
def DSquareTest_single(operator, genus, n, n_omega=11):

    deg_min = 22 - n_omega + genus - 1

    excess = 3*(genus - 1) + 2*n - 2*n_omega
    assert excess >= 0
    assert isinstance(excess, int)
    n_omega_max = n_omega + (excess // 2)
    assert isinstance(n_omega_max, int)

    GC = WOHairyGC(genus_range=[genus], 
            n_range=[n], 
            omega_range=range(n_omega, n_omega_max + 1), 
            degree_range=range(deg_min, deg_min+15), 
            differentials=[operator])

    GC.build_basis(ignore_existing_files=False)
    GC.build_matrix(ignore_existing_files=True)

    # implemented in: GraphOperator.py
    GC.square_zero_test()



def cohomology_dim_test_single(genus, n, nonzero_degrees, diagram_lists, n_omega = 11):

    assert len(nonzero_degrees) == len(diagram_lists)

    print("genus:", genus)
    print("n:", n)

    deg_min = 22 - n_omega + genus - 1

    degree_range = range(deg_min, deg_min+15)
    
    cohomdict = {}
    for degree in degree_range:

        cohom_dim = WOHairyGC.compute_cohomology_dim(degree=degree, genus=genus, n=n, n_omega=n_omega)
        
        if cohom_dim == 0:
            assert degree not in nonzero_degrees, "(g,n) = "+ str((genus, n)) + " degree: " + str(degree)
            pass
        else:
            assert cohom_dim > 0
            assert degree in nonzero_degrees, "(g,n) = "+ str((genus, n)) + " degree: " + str(degree)
            index = nonzero_degrees.index(degree)
            diagram_list = diagram_lists[index]

            test_dim = 0
            for digaram in diagram_list:
                test_dim += StandardTableaux(digaram).cardinality()
            
            print("degree:", degree, "---")
            print("actual dimension:", test_dim)
            print("nummerical dimension:", cohom_dim)
            assert cohom_dim == test_dim, "(g,n) = "+ str((genus, n))

            """
            print("d:", d)
            print("r1:", r1)
            print("r2:", r2)
            print("r3:", r3)
            print("degree:", degree)
            print("dimension:", cohom_dim)
            """

        cohomdict[degree] = cohom_dim

    print("Cohomology Dimensions (genus, n) ",
            genus, n, ":", cohomdict)



g_n_pairs = [
            (7, 2), (5, 5), (3, 8), (1, 11),
            (8, 1), (6, 4), (4, 7), (2, 10),
            (9, 0), (7, 3), (5, 6), (3, 9), (1, 12),
            (8, 2), (6, 5), (4, 8), (2, 11),
            (9, 1), (7, 4), (5, 7), (3, 10),(1, 13),
            #(10, 0), (8, 3), (9, 2), (10, 1), (11, 0)
            ]

class TestOperators(unittest.TestCase):

    def DSquareTest_EpsToOmega(self):
        for (genus, n) in g_n_pairs:
            DSquareTest_single(operator='epstoomega', genus=genus, n=n)


    def DSquareTest_ContractEdges(self):
        for (genus, n) in g_n_pairs:
            DSquareTest_single(operator='contract', genus=genus, n=n)


    def Anticommutativity_Test(self):
        for (genus, n) in g_n_pairs:
            WOHairyGC.Anticomm_Test_single(genus=genus, n=n)



    def Cohom_dim_Test(self):
        g_n_dim_pairs = [   

            # excess 0
            (7, 2, [20], [[[1,1]]]),
            (5, 5, [17], [[[1,1,1,1,1]]]),
            (3, 8, [14], [[[1,1,1,1,1,1,1,1]]]),
            (1, 11, [11], [[[1,1,1,1,1,1,1,1,1,1,1]]]),


            # excess 1
            (8, 1, [], []),
            (6, 4, [19], [[[2,1,1]]]),
            (4, 7, [16], [[[2,1,1,1,1,1]]]),
            (2, 10, [13], [[[2,1,1,1,1,1,1,1,1]]]),


            # excess 2
            (9, 0, [22], [[[]]]),
            (7, 3, [20, 21], [[[1,1,1]], [[3], [3]]]),
            (5, 6, [17, 18], [[[1,1,1,1,1,1]], [[3,1,1,1], [3,1,1,1]]]),
            (3, 9, [14, 15], [[[1,1,1,1,1,1,1,1,1]], [[3,1,1,1,1,1,1], [3,1,1,1,1,1,1]]]),
            (1, 12, [12], [[[3,1,1,1,1,1,1,1,1,1]]]),
            

            # excess 3
            (8, 2, [22], [[[2]]]),
            (6, 5, [19, 20], [[[2,1,1,1]], 
                              [[4,1], [4,1], [3,2], [3,1,1]]]),
            (4, 8, [16, 17], [[[2,1,1,1,1,1,1]], 
                              [[4,1,1,1,1], [4,1,1,1,1], [3,2,1,1,1], [3,1,1,1,1,1]]]),
            (2, 11, [14], [[[4,1,1,1,1,1,1,1], [3,2,1,1,1,1,1,1], [3,1,1,1,1,1,1,1,1]]]),


            # excess 4
            (9, 1, [22], [[[1]]]), # TODO: Mistake in Paper: "k=24 instead of k=22" ?
            (7, 4, [21, 22], [[[3,1], [3,1]], 
                              [[2,2], [2,2]]]),
            (5, 7, [18, 19], [[[3,1,1,1,1], [3,1,1,1,1]], 
                              [[5,1,1], [5,1,1], [5,1,1], [4,2,1], [4,1,1,1], [3,3,1], [3,3,1], [3,2,1,1], [3,2,1,1], [2,2,1,1,1], [2,2,1,1,1]]]),
            (3, 10, [15, 16], [[[3,1,1,1,1,1,1,1]], 
                              [[5,1,1,1,1,1], [5,1,1,1,1,1], [4,2,1,1,1,1], [4,1,1,1,1,1,1], [3,3,1,1,1,1], [3,3,1,1,1,1], [3,2,1,1,1,1,1], [3,2,1,1,1,1,1], [2,2,1,1,1,1,1,1], [2,2,1,1,1,1,1,1]]]),
            (1, 13, [13], [[[5,1,1,1,1,1,1,1,1], [3,3,1,1,1,1,1,1,1], [3,2,1,1,1,1,1,1,1,1], [2,2,1,1,1,1,1,1,1,1,1]]]),

        ]
        
        for (genus, n, nonzero_degrees, diagram_lists) in g_n_dim_pairs:
            cohomology_dim_test_single(genus=genus, n=n, nonzero_degrees=nonzero_degrees, diagram_lists=diagram_lists)
        
    
def suite():
    suite = unittest.TestSuite()
    #suite.addTest(TestOperators('DSquareTest_EpsToOmega')) # TODO: include assertion for failures
    #suite.addTest(TestOperators('DSquareTest_ContractEdges')) # TODO: include assertion for failures
    suite.addTest(TestOperators('Anticommutativity_Test'))
    suite.addTest(TestOperators('Cohom_dim_Test'))
    return suite


if __name__ == '__main__':
    print("\n#######################################\n" + "----- Start test suite WOHairy Operators -----")
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())






