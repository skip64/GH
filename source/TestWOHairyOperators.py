from WOHairyOperators import WOHairyGC

GC = WOHairyGC(genus_range=range(7,9), 
               n_range=range(0,3), 
               omega_range=range(11,13), 
               degree_range=range(10,30), 
               differentials=['epstoomega'])

# GC.build_basis(ignore_existing_files=False)
GC.build_basis(ignore_existing_files=False)

# GC.build_matrix(ignore_existing_files=False)
GC.build_matrix(ignore_existing_files=True)

GC.compute_rank(sage="integer", ignore_existing_files=False)

GC.print_dim_and_eulerchar()
GC.print_cohomology_dim()

