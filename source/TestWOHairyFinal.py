import WOHairyFinalGraphComplex
from sage.all import *

def test_basis_len(test_name, test_basis_len, genus, n, n_omega, degree):

    V = WOHairyFinalGraphComplex.WOHairyFinalGVS(genus=genus, n=n, n_omega=n_omega, degree=degree)

    if test_basis_len > 0: assert V.is_valid(), test_name

    V.build_basis(ignore_existing_files=True)

    assert V.get_dimension() == test_basis_len, test_name

    print("passed test:", test_name)


# testing excess 0 

test_basis_len("excess 0: B_1_11", 1, 1, 11, 11, 11)
test_basis_len("excess 0: B_3_8", 1, 3, 8, 11, 14)
test_basis_len("excess 0: B_5_5", 1, 5, 5, 11, 17)
test_basis_len("excess 0: B_7_2", 1, 7, 2, 11, 20)


# testing excess 1

test_basis_len("excess 1: B_2_10", 10, 2, 10, 11, 13)

test_basis_len("excess 1: B_6_4", 5, 6, 4, 11, 19)


# testing excess 2

test_basis_len("excess 2: B_7_3", 16, 7, 3, 11, 20)


# testing excess 3

V = WOHairyFinalGraphComplex.WOHairyFinalGVS(genus=8, n=2, n_omega=11, degree=20)

print("is valid: " + str(V.is_valid()))
print("excess:", V.excess)

# Compute basis
V.build_basis(ignore_existing_files=False)

print("dimension:", V.get_dimension())

# Plot the basis elements to images, open html page (temp/temp.html) with the images
print("now plotting ----")
V.display_basis_plots()