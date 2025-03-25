from sage.all import *

G = Graph({0:[3],1:[3],2:[3]})

print("before relabelling: ---")
print(G.vertices())
print(G.edges())


#G = G.relabel({0:5, 1:4, 2:3, 3:2, 4:1, 5:0}, inplace=False)
#G.relabel([5,4,3,2,1,0], inplace=True)
G.relabel([1,2,0], inplace=True)

print("after relabelling: ---")
print(G.vertices())
#print(G.vertices(sort=True))
print(G.edges())

