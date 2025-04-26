# testing of operator-functionality
# for testing to be fully rigurous, the "data"-folder should be deleted before running the tests


import unittest
from WOHairyOperators import WOHairyGC
from sage.all import *



class TestEpsToOmega(unittest.TestCase):

    @staticmethod
    def DSquareTest_single(genus, n, n_omega, min_degree, max_degree):
        # min_degree: minimal degree with nonzero dimension
        # max_degree: maximal degree with nonzero dimension

        excess = 3*(genus - 1) + 2*n - 2*n_omega
        assert excess >= 0
        assert isinstance(excess, int)
        n_omega_max = n_omega + (excess // 2)
        assert isinstance(n_omega_max, int)

        GC = WOHairyGC(genus_range=[genus], 
               n_range=[n], 
               omega_range=range(n_omega, n_omega_max + 1), 
               degree_range=range(min_degree-1, max_degree+1), 
               differentials=['epstoomega'])

        GC.build_basis(ignore_existing_files=False)
        GC.build_matrix(ignore_existing_files=True)

        # TODO: make fail if not "success" -> Operator-Class
        GC.square_zero_test()


    # TODO: add more test-cases
    def DSquareTest(self):
        self.DSquareTest_single(genus=9, n=1, n_omega=11, min_degree=21, max_degree=25)
        self.DSquareTest_single(genus=8, n=3, n_omega=11, min_degree=19, max_degree=24)
        self.DSquareTest_single(genus=11, n=0, n_omega=11, min_degree=22, max_degree=30)

        
    
def suite():
    suite = unittest.TestSuite()
    suite.addTest(TestEpsToOmega('DSquareTest'))
    return suite


if __name__ == '__main__':
    print("\n#######################################\n" + "----- Start test suite WOHairy Operators -----")
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())






