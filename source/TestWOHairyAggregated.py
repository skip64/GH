from WOHairyGC import WOHairyAggregatedGVS



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


V = WOHairyAggregatedGVS(n_components=n_components, 
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