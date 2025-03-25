import WOHairyFinalGraphComplex
from sage.all import *


def compute_euler_char(genus, n):

    euler_char = 0

    for n_omega in [11, 12, 13, 14, 15, 16]:
        
        omega_excess = 3*(genus-1) + 2*n - 2*n_omega

        if omega_excess < 0: break
        print("---")
        print("n_omega:", n_omega)
        print("excess:", omega_excess)

        # degree lower-bound
        # n_vertices = degree - 22 + n_omega - genus + 1 >= 0
        # -> degree >= 22 - n_omega + genus - 1
        deg_min = 22 - n_omega + genus - 1

        euler_char_omega = 0

        for degree in range(deg_min, deg_min+15):
            
            #print(genus, n, n_omega, degree)
            V = WOHairyFinalGraphComplex.WOHairyFinalGVS(genus=genus, n=n, n_omega=n_omega, degree=degree)

            V.build_basis(ignore_existing_files=False)

            print("degree:", degree, " dimension:", V.get_dimension())
            
            euler_char_omega += (-1)**degree * V.get_dimension()

        print("contribution:", euler_char_omega)
        euler_char += euler_char_omega

    return euler_char


def test_euler_char(genus, n, coefficients, diagrams, excess):

    print("testing for (g,n)=", (genus, n), "-----------------------")

    computed_excess = 3*genus + 2*n - 25

    if computed_excess < 0: assert excess < 0
    else: assert computed_excess == excess

    euler_char = compute_euler_char(genus, n)
    
    assert len(coefficients) == len(diagrams)

    euler_char_reference = 0
    for coefficient, diagram in zip(coefficients, diagrams):
        euler_char_reference += coefficient * StandardTableaux(diagram).cardinality()

    print("euler_char:", euler_char)
    print("euler_char_reference:", euler_char_reference)
    assert euler_char == euler_char_reference, (euler_char, euler_char_reference)


# excess < 0
for i in range(5): test_euler_char(5, i, coefficients=[], diagrams=[], excess=-1)
for i in range(4): test_euler_char(6, i, coefficients=[], diagrams=[], excess=-1)
for i in range(2): test_euler_char(7, i, coefficients=[], diagrams=[], excess=-1)
test_euler_char(8, 0, coefficients=[], diagrams=[], excess=-1)

# excess = 0
test_euler_char(5, 5, coefficients=[-1], diagrams=[[1,1,1,1,1]], excess=0)
test_euler_char(7, 2, coefficients=[1], diagrams=[[1,1]], excess=0)

# excess = 1
test_euler_char(6, 4, coefficients=[-1], diagrams=[[2,1,1]], excess=1)
test_euler_char(8, 1, coefficients=[], diagrams=[], excess=1)

# excess = 2
test_euler_char(5, 6, coefficients=[-1, 2], diagrams=[[1,1,1,1,1,1], [3,1,1,1]], excess=2)
test_euler_char(7, 3, coefficients=[1, -2], diagrams=[[1,1,1], [3]], excess=2)
test_euler_char(9, 0, coefficients=[1], diagrams=[[]], excess=2)

# excess = 3
test_euler_char(8, 2, coefficients=[1], diagrams=[[2]], excess=3)
test_euler_char(6, 5, coefficients=[-1, 1, 1, 2], diagrams=[[2,1,1,1], [3,1,1], [3,2], [4,1]], excess=3)

# excess = 4
test_euler_char(9, 1, coefficients=[1], diagrams=[[1]], excess=4)
#test_euler_char(7, 4, coefficients=[2, -2], diagrams=[[2,2], [3,1]], excess=4)

# excess = 5
test_euler_char(10, 0, coefficients=[-2], diagrams=[[]], excess=5)
#test_euler_char(8, 3, coefficients=[4], diagrams=[[1,1,1]], excess=5)
#test_euler_char(6, 6, coefficients=[-4, -2, -3, -3, 1, -2, -2, 2, -3, -2 ,-3], 
#                diagrams=[[1,1,1,1,1,1], [2,1,1,1,1], [2,2,1,1], [2,2,2], 
#                          [3,1,1,1], [3,2,1], [3,3], [4,1,1], [4,2], [5,1], [6]], excess=5)

# excess = 6
#test_euler_char(9, 2, coefficients=[4], diagrams=[[2]], excess=6)
#test_euler_char(7, 5, coefficients=[-1, -7, -6, 2, 1, 4], diagrams=[[1,1,1,1,1], [2,1,1,1], [3,1,1], [3,2], [4,1], [5]], excess=6)

# excess = 7
#test_euler_char(10, 1, coefficients=[-6], diagrams=[[1]], excess=7)
#test_euler_char(8, 4, coefficients=[10, 2, -2, -10, -3], diagrams=[[1,1,1,1], [2,1,1], [2,2], [3,1], [4]], excess=7)


# excess = 8
#test_euler_char(11, 0, coefficients=[2], diagrams=[[]], excess=8)
#test_euler_char(9, 3, coefficients=[1, 11, 4], diagrams=[[1,1,1], [2,1], [3]], excess=8)

# excess = 9
#test_euler_char(10, 2, coefficients=[-9, -2], diagrams=[[1,1], [2]], excess=9)

# excess = 10
#test_euler_char(11, 1, coefficients=[1], diagrams=[[1]], excess=10)

