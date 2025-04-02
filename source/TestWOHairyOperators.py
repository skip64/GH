from WOHairyOperators import WOHairyGC

# Testing H^d(5,6) i.e. genus = 5, n = 6
# nontrivial for d=17,18
# excess = 2 so omega = 11,12
GC = WOHairyGC(genus_range=range(5,6), 
               n_range=range(6,7), 
               omega_range=range(11,13), 
               degree_range=range(15,23), 
               differentials=['epstoomega', 'contract'])


GC.build_basis(ignore_existing_files=False)

GC.build_matrix(ignore_existing_files=True)

GC.compute_rank(sage="integer", ignore_existing_files=True)

GC.print_dim_and_eulerchar()
GC.print_cohomology_dim()

