from WOHairyGC import WOHairyGC
import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.preprocessing import PolynomialFeatures
from sklearn.model_selection import cross_val_score


"""
g_n_pairs = [
    (7, 2), (5, 5), (3, 8), (1, 11),
    (8, 1), (6, 4), (4, 7), (2, 10), (0, 13),
    (9, 0), (7, 3), (5, 6), (3, 9), (1, 12),
    (8, 2), (6, 5), (4, 8), (2, 11), (0, 14),
    (9, 1), (7, 4), (5, 7), (3, 10), (1, 13),
    (10, 0), (8, 3), (6, 6),
    (9, 2), (7, 5),
    (10, 1), 
    (11, 0)
    ]
"""


g_n_pairs = [(0, 0), (0, 1), (1, 0), (0, 2), (1, 1), (2, 0), (0, 3), (1, 2), (2, 1), (3, 0), (0, 4), (1, 3), (2, 2), (3, 1), (4, 0), (0, 5), (1, 4), (2, 3), (3, 2), (4, 1), (0, 6), (5, 0), (1, 5), (2, 4), (3, 3), (4, 2), (0, 7), (5, 1), (1, 6), (6, 0), (2, 5), (3, 4), (4, 3), (0, 8), (5, 2), (1, 7), (6, 1), (2, 6), (7, 0), (3, 5), (4, 4), (0, 9), (5, 3), (1, 8), (6, 2), (2, 7), (7, 1), (3, 6), (8, 0), (4, 5), (0, 10), (5, 4), (1, 9), (6, 3), (2, 8), (7, 2), (3, 7), (8, 1), (4, 6), (0, 11), (9, 0), (5, 5), (1, 10), (6, 4), (2, 9), (7, 3), (3, 8), (8, 2), (4, 7), (0, 12), (9, 1), (5, 6), (1, 11), (10, 0), (6, 5), (2, 10), (7, 4), (3, 9), (8, 3), (4, 8), (0, 13), (9, 2), (5, 7), (1, 12), (10, 1), (6, 6), (2, 11), (11, 0), (7, 5), (3, 10), (8, 4), (4, 9), (0, 14), 
             #(9, 3), (5, 8), (1, 13), (10, 2), (6, 7), (2, 12), (11, 1), (7, 6), (3, 11), (12, 0)
             ]


max_dimensions = []
relevant_g_n_pairs = []

for (genus, n) in g_n_pairs:
    max_dimension = WOHairyGC.max_basis_dimension(genus=genus, n=n)
    if max_dimension > 100:
        max_dimensions.append(max_dimension)
        relevant_g_n_pairs.append((genus, n))

log_max_dimensions = np.log(max_dimensions)

print(log_max_dimensions)

X = np.array(relevant_g_n_pairs, dtype=float)
y = log_max_dimensions

model = LinearRegression()
model.fit(X, y)

print("Coefficients:", model.coef_)
print("Intercept:", model.intercept_)
print("R^2 score:", model.score(X, y))


def predict(genus, n):
    return np.exp(model.predict(np.array([[genus, n]]))[0])



for i, (genus, n) in enumerate(relevant_g_n_pairs):
    print(f"Predicted dimension for genus {genus}, n {n}: {predict(genus, n)}")
    #print(WOHairyGC.max_basis_dimension_estimate(excess, n))
    print(f"Actual dimension for genus {genus}, n {n}: {max_dimensions[i]}")
    






