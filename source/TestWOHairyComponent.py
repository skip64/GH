from WOHairyBasisGeneration import WOHairyComponentGVS


# Instance with 5 vertices, 2 loops, 4 hairs
#V = _ModifiedCHairyGraphComplex.ModifiedCHairyGVS(n_vertices=5, n_loops=2, n=1, n_omega=2, n_epsilon=1)

#V = _ModifiedCHairyGraphComplex.ModifiedCHairyGVS(n_vertices=3, n_loops=1, n=1, n_omega=2, n_epsilon=1)

#V = _ModifiedCHairyGraphComplex.ModifiedCHairyGVS(n_vertices=1, n_loops=0, n=1, n_omega=1, n_epsilon=1)

# For testing perm_sign bad behaviour!
#V = _ModifiedCHairyGraphComplex.ModifiedCHairyGVS(n_vertices=1, n_loops=0, n=0, n_omega=0, n_epsilon=3)

#V = _ModifiedCHairyGraphComplex.ModifiedCHairyGVS(n_vertices=0, n_loops=0, n=1, n_omega=0, n_epsilon=1)

# 
V = WOHairyComponentGVS(n_vertices=2, n_loops=0, n=2, n_omega=2, n_epsilon=0)

print("is valid: " + str(V.is_valid()))

#print(V.get_ordered_param_dict().get_value_tuple())

# Compute basis
V.build_basis(ignore_existing_files=True)

# Iterate over basis, displaying every vector (as g6 ascii)
#for s in V.get_basis_g6():
#      print(s)

# Plot the basis elements to images, open html page (temp/temp.html) with the images
V.display_basis_plots()
