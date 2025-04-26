from WOHairyOperators import WOHairyGC



# Testing H^d(9,1) i.e. genus = 9, n = 1
# nontrivial for d = 21,22,23,24,25
# excess = 4 so omega = 11,12,13
GC = WOHairyGC(genus_range=range(9,10), 
               n_range=range(1,2), 
               omega_range=range(11,14), 
               degree_range=range(20,26), 
               differentials=['epstoomega', 'contract'])


GC.build_basis(ignore_existing_files=False)

GC.build_matrix(ignore_existing_files=True)

GC.square_zero_test()


"""
EpsToOmegaGO ---
def is_match(domain, target):
    return (domain.genus == target.genus 
            and domain.n == target.n
            and domain.n_omega + 1 == target.n_omega 
            and domain.degree - 1 == target.degree)

ContractEdgesGO ---
def is_match(domain, target):
    return (domain.genus == target.genus 
            and domain.n == target.n
            and domain.n_omega == target.n_omega 
            and domain.degree - 1 == target.degree)
"""

"""
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

GC.square_zero_test()
"""

"""
Square zero test for <contract edges differential>> failed for the pair:
<contract edges graph operator, domain: <wohairy vector space with parameters: (genus: 5 n: 6 omegas: 11 degree: 18 )>>, 
<contract edges graph operator, domain: <wohairy vector space with parameters: (genus: 5 n: 6 omegas: 11 degree: 17 )>>
"""

#GC.compute_rank(sage="integer", ignore_existing_files=True)

#GC.print_dim_and_eulerchar()
#GC.print_cohomology_dim()






