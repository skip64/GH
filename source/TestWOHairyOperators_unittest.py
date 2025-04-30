# testing of operator-functionality
# for testing to be fully rigurous, the "data"-folder should be deleted before running the tests


import unittest
from WOHairyOperators import WOHairyGC
from sage.all import *
import GraphOperator
import itertools
import Parameters
import Shared

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


def Anticomm_Test(genus, n, n_omega=11, eps=Parameters.square_zero_test_eps):

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
            differentials=['contract', 'epstoomega'])

    GC.build_basis(ignore_existing_files=False)
    GC.build_matrix(ignore_existing_files=True)

    print(GC.operator_collection_list)
    assert len(GC.operator_collection_list) == 2

    dif_1 = GC.operator_collection_list[0]
    dif_2 = GC.operator_collection_list[1]
    op1_matrix_list = dif_1.op_matrix_list
    op2_matrix_list = dif_2.op_matrix_list
    
    for (op1_1, op1_2, op2_1, op2_2) in itertools.product(op1_matrix_list, op1_matrix_list, op2_matrix_list, op2_matrix_list):
        if op1_1.get_domain() == op2_2.get_domain() and op1_2.get_target() == op2_1.get_target() \
            and op1_1.get_target() == op2_1.get_domain() and op2_2.get_target() == op1_2.get_domain():
            if (not (op1_1.is_valid() and op1_2.is_valid() and op2_1.is_valid() and op2_2.is_valid())) \
                or ((op1_1.is_trivial() or op2_1.is_trivial()) and (op1_2.is_trivial() or op2_2.is_trivial())):
                #print('trivial success')
                pass

            else: 
                M1_1 = op1_1.get_matrix()
                M1_2 = op1_2.get_matrix()
                M2_1 = op2_1.get_matrix()
                M2_2 = op2_2.get_matrix()

                product_1 = M2_1 * M1_1
                product_2 = M1_2 * M2_2

                if Shared.matrix_norm(product_1 - product_2) < eps:
                    print('success')
                else:
                    assert False, 'Anitcomm-Test failed!'


g_n_pairs = [(5, 5), (7, 2), (6, 4), (8, 1), 
                (5, 6), (7, 3), (9, 0), (8, 2), (6, 5),
                (9, 1), (7, 4), (10, 0), (8, 3),
                (9, 2), (10, 1), (11, 0)]

class TestOperators(unittest.TestCase):

    def DSquareTest_EpsToOmega(self):
        for (genus, n) in g_n_pairs:
            DSquareTest_single(operator='epstoomega', genus=genus, n=n)


    def DSquareTest_ContractEdges(self):
        for (genus, n) in g_n_pairs:
            DSquareTest_single(operator='contract', genus=genus, n=n)


    def Anticommutativity_Test(self):
        for (genus, n) in g_n_pairs:
            Anticomm_Test(genus=genus, n=n)
        
    
def suite():
    suite = unittest.TestSuite()
    #suite.addTest(TestOperators('DSquareTest_EpsToOmega'))
    #suite.addTest(TestOperators('DSquareTest_ContractEdges'))
    suite.addTest(TestOperators('Anticommutativity_Test'))
    return suite


if __name__ == '__main__':
    print("\n#######################################\n" + "----- Start test suite WOHairy Operators -----")
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())






