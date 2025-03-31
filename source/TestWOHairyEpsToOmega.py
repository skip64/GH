from WOHairyOperators import EpsToOmegaGO

genus = 9
n = 0
n_omega = 11
degree = 22

Operator = EpsToOmegaGO.generate_operator(genus, n, n_omega, degree)

D = Operator.get_matrix()

