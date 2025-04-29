# testing of operator-functionality
# for testing to be fully rigurous, the "data"-folder should be deleted before running the tests


import unittest
from WOHairyOperators import WOHairyGC
from sage.all import *

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

    # TODO: make fail if not "success" -> Operator-Class
    GC.square_zero_test()


class TestOperators(unittest.TestCase):

    def DSquareTest_EpsToOmega(self):

        g_n_pairs = [(5, 5), (7, 2), (6, 4), (8, 1), 
                     (5, 6), (7, 3), (9, 0), (8, 2), (6, 5),
                     (9, 1), (7, 4), (10, 0), (8, 3),
                     (9, 2), (10, 1), (11, 0)]
        
        for (genus, n) in g_n_pairs:
            DSquareTest_single(operator='epstoomega', genus=genus, n=n)


    def DSquareTest_ContractEdges(self):

        g_n_pairs = [(5, 5), (7, 2), (6, 4), (8, 1), 
                     (5, 6), (7, 3), (9, 0), (8, 2), (6, 5),
                     (9, 1), (7, 4), (10, 0), (8, 3),
                     (9, 2), (10, 1), (11, 0)]
        
        for (genus, n) in g_n_pairs:
            DSquareTest_single(operator='contract', genus=genus, n=n)

        
    
def suite():
    suite = unittest.TestSuite()
    #suite.addTest(TestOperators('DSquareTest_EpsToOmega'))
    suite.addTest(TestOperators('DSquareTest_ContractEdges'))
    return suite


if __name__ == '__main__':
    print("\n#######################################\n" + "----- Start test suite WOHairy Operators -----")
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())






