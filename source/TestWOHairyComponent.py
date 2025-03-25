import WOHairyComponentGraphComplex
from sage.all import *


def test_basis_len(test_name, test_basis_len, n_vertices, n_loops, n, n_omega, n_epsilon):

    V = WOHairyComponentGraphComplex.WOHairyComponentGVS(n_vertices=n_vertices, n_loops=n_loops, n=n, n_omega=n_omega, n_epsilon=n_epsilon)

    if test_basis_len > 0: assert V.is_valid(), test_name

    V.build_basis(ignore_existing_files=True)

    actual_basis_len = len([G for G in V.get_basis()])
    assert actual_basis_len == test_basis_len, test_name
    
    print("passed test:", test_name)


# Testing double-leg configurations

test_basis_len("omega-omega", test_basis_len=0, n_vertices=0, n_loops=0, n=0, n_omega=2, n_epsilon=0)

test_basis_len("epsilon-epsilon", test_basis_len=1, n_vertices=0, n_loops=0, n=0, n_omega=0, n_epsilon=2)

test_basis_len("omega-epsilon", test_basis_len=1, n_vertices=0, n_loops=0, n=0, n_omega=1, n_epsilon=1)

test_basis_len("leg-omega", test_basis_len=1, n_vertices=0, n_loops=0, n=1, n_omega=1, n_epsilon=0)

test_basis_len("leg-epsilon", test_basis_len=1, n_vertices=0, n_loops=0, n=1, n_omega=0, n_epsilon=1)


# testing for symmetries in graphs with one vertex

test_basis_len("eps-omega-eps", test_basis_len=0, n_vertices=1, n_loops=0, n=0, n_omega=1, n_epsilon=2)

test_basis_len("3-omega", test_basis_len=1, n_vertices=1, n_loops=0, n=0, n_omega=3, n_epsilon=0)

test_basis_len("3-epsilon", test_basis_len=0, n_vertices=1, n_loops=0, n=0, n_omega=0, n_epsilon=3)




# Instance with 5 vertices, 2 loops, 4 hairs
#V = _ModifiedCHairyGraphComplex.ModifiedCHairyGVS(n_vertices=5, n_loops=2, n=1, n_omega=2, n_epsilon=1)

#V = _ModifiedCHairyGraphComplex.ModifiedCHairyGVS(n_vertices=3, n_loops=1, n=1, n_omega=2, n_epsilon=1)

#V = _ModifiedCHairyGraphComplex.ModifiedCHairyGVS(n_vertices=1, n_loops=0, n=1, n_omega=1, n_epsilon=1)

# For testing perm_sign bad behaviour!
#V = _ModifiedCHairyGraphComplex.ModifiedCHairyGVS(n_vertices=1, n_loops=0, n=0, n_omega=0, n_epsilon=3)

#V = _ModifiedCHairyGraphComplex.ModifiedCHairyGVS(n_vertices=0, n_loops=0, n=1, n_omega=0, n_epsilon=1)

# 
V = WOHairyComponentGraphComplex.WOHairyComponentGVS(n_vertices=2, n_loops=0, n=2, n_omega=2, n_epsilon=0)

print("is valid: " + str(V.is_valid()))

#print(V.get_ordered_param_dict().get_value_tuple())

# Compute basis
V.build_basis(ignore_existing_files=True)

# Iterate over basis, displaying every vector (as g6 ascii)
#for s in V.get_basis_g6():
#      print(s)

# Plot the basis elements to images, open html page (temp/temp.html) with the images
V.display_basis_plots()
