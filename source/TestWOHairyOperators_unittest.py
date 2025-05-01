# testing of operator-functionality
# for testing to be fully rigurous, the "data"-folder should be deleted before running the tests


import unittest
from WOHairyOperators import WOHairyGC, ContractEdgesGO, EpsToOmegaGO, WOHairyFinalGVS
from sage.all import *
import GraphOperator
import itertools
import Parameters
import Shared
from sage.all import StandardTableaux


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


def Anticomm_Test_single(genus, n, n_omega=11, eps=Parameters.square_zero_test_eps):

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



def cohomology_dim_test_single(genus, n, degrees, diagrams, n_omega = 11):

    print("genus:", genus)
    print("n:", n)

    deg_min = 22 - n_omega + genus - 1

    degree_range = range(deg_min, deg_min+15)
    
    cohomdict = {}
    for degree in degree_range:

        D_Cont_deg = ContractEdgesGO.generate_operator(degree=degree, genus=genus, n=n, n_omega=n_omega)
        D_Cont_deg_p1 = ContractEdgesGO.generate_operator(degree=degree+1, genus=genus, n=n, n_omega=n_omega)
        D_eps_deg = EpsToOmegaGO.generate_operator(degree=degree, genus=genus, n=n, n_omega=n_omega)
        D_eps_deg_p1 = EpsToOmegaGO.generate_operator(degree=degree+1, genus=genus, n=n, n_omega=n_omega)

        for op in [D_Cont_deg, D_Cont_deg_p1, D_eps_deg, D_eps_deg_p1]:
            op.domain.build_basis()
            op.target.build_basis()

            op.build_matrix(skip_if_no_basis=False)

        D_Cont_deg_mat = D_Cont_deg.get_matrix()
        D_Cont_deg_p1_mat = D_Cont_deg_p1.get_matrix()
        D_eps_deg_mat = D_eps_deg.get_matrix()
        D_eps_deg_p1_mat = D_eps_deg_p1.get_matrix()

        assert D_Cont_deg.domain == D_eps_deg.domain
        deg_double_mat = D_Cont_deg_mat.stack(D_eps_deg_mat)

        assert D_Cont_deg_p1.domain == D_eps_deg_p1.domain
        deg_p1_double_mat = D_Cont_deg_p1_mat.stack(D_eps_deg_p1_mat) 
        
        d = WOHairyFinalGVS(genus=genus, n=n, n_omega=n_omega, degree=degree).get_dimension()
        r1 = deg_double_mat.rank()
        r2 = deg_p1_double_mat.rank()
        r3 = D_eps_deg_p1_mat.rank()

        cohom_dim = d - r1 - r2 + r3
        
        if cohom_dim > 0:
            print("d:", d)
            print("r1:", r1)
            print("r2:", r2)
            print("r3:", r3)
            print("degree:", degree)
            print("dimension:", d - r1 - r2 + r3)

            cohomdict[degree] = cohom_dim

    print("Cohomology Dimensions (genus, n) ",
            genus, n, ":", cohomdict)
    
    for degree, diagram in zip(degrees, diagrams):
        assert cohomdict[degree] == StandardTableaux(diagram).cardinality()


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
            Anticomm_Test_single(genus=genus, n=n)



    def Cohom_dim_Test(self):
        g_n_dim_pairs = [   # excess 0
                            (1, 11, [11], [[1,1,1,1,1,1,1,1,1,1,1]]),
                            (3, 8, [14], [[1,1,1,1,1,1,1,1]]),
                            (5, 5, [17], [[1,1,1,1,1]]),
                            (7, 2, [20], [[1,1]]),
                            # excess 1
                            (2, 10, [13], [[2,1,1,1,1,1,1,1,1]]),
                            (4, 7, [16], [[2,1,1,1,1,1]])
                         ]
        for (genus, n, degrees, diagrams) in g_n_dim_pairs:
            cohomology_dim_test_single(genus=genus, n=n, degrees=degrees, diagrams=diagrams)
        
    
def suite():
    suite = unittest.TestSuite()
    #suite.addTest(TestOperators('DSquareTest_EpsToOmega'))
    #suite.addTest(TestOperators('DSquareTest_ContractEdges'))
    #suite.addTest(TestOperators('Anticommutativity_Test'))
    suite.addTest(TestOperators('Cohom_dim_Test'))
    return suite


if __name__ == '__main__':
    print("\n#######################################\n" + "----- Start test suite WOHairy Operators -----")
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())






