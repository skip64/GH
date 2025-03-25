import WOHairyAggregatedGraphComplex
from sage.all import *
import math


def test_basis_len(test_name, test_basis_len, inner_edges, n_components, n_vertices, genus, n, n_omega, n_epsilon, n_double_legs):

    V = WOHairyAggregatedGraphComplex.WOHairyAggregatedGVS(n_components=n_components, n_vertices=n_vertices, genus=genus, n=n, n_omega=n_omega, n_epsilon=n_epsilon, n_double_legs=n_double_legs, do_print=False)

    assert V.n_edges == inner_edges, inner_edges

    if test_basis_len > 0: assert V.is_valid(), test_name

    V.build_basis(ignore_existing_files=True)

    actual_basis_len = len([G for G in V.get_basis()])
    assert actual_basis_len == test_basis_len, test_name

    print("passed test:", test_name)



# Testing double-leg configurations ---


n_components = 1
n_vertices = 0
n_double_legs = 1

test_basis_len("omega-omega", test_basis_len=0, inner_edges=0, genus=2, n=0, n_omega=2, n_epsilon=0, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)

test_basis_len("epsilon-epsilon", test_basis_len=1, inner_edges=0, genus=2, n=0, n_omega=0, n_epsilon=2, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)

test_basis_len("omega-epsilon", test_basis_len=1, inner_edges=0, genus=2, n=0, n_omega=1, n_epsilon=1, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)

test_basis_len("leg-omega", test_basis_len=1, inner_edges=0, genus=1, n=1, n_omega=1, n_epsilon=0, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)

test_basis_len("leg-epsilon", test_basis_len=1, inner_edges=0, genus=1, n=1, n_omega=0, n_epsilon=1, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)




# testing connected graphs with one vertex ---


n_components = 1
n_vertices = 1
n_double_legs = 0

test_basis_len("eps-omega-eps", test_basis_len=0, inner_edges=0, genus=3, n=0, n_omega=1, n_epsilon=2, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)

test_basis_len("3-omega", test_basis_len=1, inner_edges=0, genus=3, n=0, n_omega=3, n_epsilon=0, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)

test_basis_len("3-epsilon", test_basis_len=0, inner_edges=0, genus=3, n=0, n_omega=0, n_epsilon=3, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs)




# testing graphs with multiple components ---

test_basis_len("omega-1 & omega-2", test_basis_len=1, inner_edges=0, n_components=2, n_vertices=0, genus=1, n=2, n_omega=2, n_epsilon=0, n_double_legs=2)

test_basis_len("omega-1 & omega-2 & 2x3-omega", test_basis_len=1, inner_edges=0, n_components=4, n_vertices=2, genus=5, n=2, n_omega=8, n_epsilon=0, n_double_legs=2)


# excess 0 examples ---

test_basis_len("excess 0: B_(7,2)", test_basis_len=1, inner_edges=0, n_components=5, 
               n_vertices=3, genus=7, n=2, n_omega=11, n_epsilon=0, n_double_legs=2)



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

                V = WOHairyAggregatedGraphComplex.WOHairyAggregatedGVS(n_components=n_components, 
                                                n_vertices=n_vertices, 
                                                genus=genus, 
                                                n=n, n_omega=n_omega, n_epsilon=n_epsilon,
                                                n_double_legs=n_double_legs,
                                                do_print=False)
                
                assert V.excess == excess

                V.build_basis(ignore_existing_files=True)

                #V.display_basis_plots()

                assert V.get_dimension() == 1, (excess, n_omega, n)



"""
# two epsilon double-legs
n_components = 2
n_vertices = 0
genus = 3
n = 0
n_omega = 0
n_epsilon = 4
n_double_legs = 2
"""

"""
# three epsilon double-legs
n_components = 3
n_vertices = 0
genus = 4
n = 0
n_omega = 0
n_epsilon = 6
n_double_legs = 3
"""

# rebuilding the example from the paper

"""
# eps-1 double-leg
n_components = 1
n_vertices = 0
genus = 1
n = 1
n_omega = 0
n_epsilon = 1
n_double_legs = 1
"""

"""
# non-trivial example
n_components = 1
n_vertices = 3
genus = 3
n = 1
n_omega = 1
n_epsilon = 1
n_double_legs = 0
"""

"""
# non-trivial example
n_components = 1
n_vertices = 4
genus = 4
n = 1
n_omega = 2
n_epsilon = 1
n_double_legs = 0
"""

"""
# very non-trivial example
n_components = 2
n_vertices = 7
genus = 6
n = 2
n_omega = 3
n_epsilon = 2
n_double_legs = 0
"""

"""
# example from the paper
n_components = 5
n_vertices = 4
genus = 9
n = 1
n_omega = 12
n_epsilon = 1
n_double_legs = 1
"""


# double 3-omega
n_components = 2
n_vertices = 2
genus = 5
n = 0
n_omega = 6
n_epsilon = 0
n_double_legs = 0


# double 3-omega with omega-1 leg
n_components = 3
n_vertices = 2
genus = 5
n = 1
n_omega = 7
n_epsilon = 0
n_double_legs = 1

# two omega-j legs
n_components = 2
n_vertices = 0
genus = 1
n = 2
n_omega = 2
n_epsilon = 0
n_double_legs = 2

"""
"""

# very non-trivial example
n_components = 2
n_vertices = 7
genus = 6
n = 2
n_omega = 3
n_epsilon = 2
n_double_legs = 0


# excess 0 example in B_(11,1)
n_components = 11
n_vertices = 0
genus = 1
n = 11
n_omega = 11
n_epsilon = 0
n_double_legs = 11


#problematic graph:
n_components = 7
n_vertices = 2
genus=6
n=5
n_omega=11
n_epsilon = 0
n_double_legs = 5

#problematic graph 2:
n_components = 6
n_vertices = 2
genus=6
n=4
n_omega=10
n_epsilon = 0
n_double_legs = 4

#problematic graph 3:
n_components = 5
n_vertices = 2
genus=6
n=3
n_omega=9
n_epsilon = 0
n_double_legs = 3

#problematic graph 4:
n_components = 4
n_vertices = 2
genus=6
n=2
n_omega=8
n_epsilon = 0
n_double_legs = 2

n_components = 3
n_vertices = 2
genus=6
n=1
n_omega=7
n_epsilon = 0
n_double_legs = 1


n_components = 2
n_vertices = 2
genus=6
n=0
n_omega=6
n_epsilon = 0
n_double_legs = 0


#problematic graph:
n_components = 7
n_vertices = 2
genus=6
n=5
n_omega=11
n_epsilon = 0
n_double_legs = 5

#problematic graph:
n_components = 6
n_vertices = 2
genus=6
n=4
n_omega=10
n_epsilon = 0
n_double_legs = 4

n_components = 5
n_vertices = 2
genus=6
n=3
n_omega=9
n_epsilon = 0
n_double_legs = 3

n_components = 3
n_vertices = 2
genus=6
n=1
n_omega=7
n_epsilon = 0
n_double_legs = 1





                



n_components = 1
n_vertices = 2
genus=2
n=2
n_omega=2
n_epsilon = 0
n_double_legs = 0


V = WOHairyAggregatedGraphComplex.WOHairyAggregatedGVS(n_components=n_components, 
                                                n_vertices=n_vertices, 
                                                genus=genus, 
                                                n=n, n_omega=n_omega, n_epsilon=n_epsilon,
                                                n_double_legs=n_double_legs,
                                                do_print=True)

print("excess", V.excess)

#print(V.get_ordered_param_dict().get_value_tuple())

# Compute basis
V.build_basis(ignore_existing_files=True)

# Iterate over basis, displaying every vector (as g6 ascii)
#for s in V.get_basis_g6():
#      print(s)

# Plot the basis elements to images, open html page (temp/temp.html) with the images
V.display_basis_plots()