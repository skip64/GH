from WOHairyGC_Pascal import WOHairyFinalGVS

# Manual Testing

# big components
#V = WOHairyFinalGraphComplex.WOHairyFinalGVS(genus=11, n=0, n_omega=11, degree=30)

V = WOHairyFinalGVS(genus=9, n=0, n_omega=11, degree=22)

print("is valid: " + str(V.is_valid()))
print("excess:", V.excess)

# Compute basis
V.build_basis(ignore_existing_files=False)

print("dimension:", V.get_dimension())

# Plot the basis elements to images, open html page (temp/temp.html) with the images
print("now plotting ----")
V.display_basis_plots()