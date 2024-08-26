# ---------------------------------------------------------------- 
# Written by: Hyoung Suk Suh (h.suh@columbia.edu)     
# sub by: Chulmin Kweon      (Chulmin.Kweon@columbia.edu)
# ----------------------------------------------------------------       

import numpy as np
import sys
sys.path.append("./include")
import time as ct

import Program
import WriteResult
import ElastoPlasticity as Plasticity
import FemBulk as Bulk
from MeshHandle import *
import importlib



max_itr = 20 # Global max iteration
tic = ct.time()
assert len(sys.argv) == 2, "check the input value\n execution example:\n\t python3 MainElastoPlastic.py ./input/PlateWithHole/PlateWithHole_113_93_plastic.inp"


# Model Generation
Fem = Program.Model()
input_name = sys.argv[1]

# Read input
Node, Element = Program.ReadInput(input_name, Fem)
#print ("FEM analysis is starting in 2 seconds")
#ct.sleep(2)
EdgeNode    = []
EdgeElement = []

##print(Node.NNode)
#Element.Connectivity = [[1,2,3,4],
                        #[2,5,6,3]]
#print(Element.Connectivity)
#exit(1)
#for kk, elem in enumerate(Element.Connectivity):
    #elem = elem.copy()
    #elem.append(elem[0])
    #for ii in range(len(elem)-1):
        ##pri
        #flag = 0
        #for jj, EN in enumerate(EdgeNode):
            #if EN[0] == elem[ii+1] and EN[1] == elem[ii]:
                #flag = jj
                #break

        #if flag == 0:
            #EdgeNode.append([elem[ii], elem[ii+1]])
            #EdgeElement.append([kk,-1])
        #else:
            #print(jj)
            #print(EdgeElement)
            #EdgeElement[jj][1] = kk
            #print(EdgeElement)
            #print(kk)
            #input("")

#print(EdgeNode)
#print(len(EdgeNode))
#print("")
#print(EdgeElement)
#print(EdgeElement[1][1])
#print(len(EdgeElement))
#print("")
#print(Element.Connectivity)
#print("")
Element.Connectivity = [[1,2,3,4],
                        [2,5,6,3],
                        [2,7,8,5],
                        [9,7,2,1]]
Element.NElem = 4
Element.Id = [1,2,3,4]
Edge = GenerateEdge(Element)
for ii in range(len(Edge.AdjacElem)):
    print(Edge.AdjacElem[ii], "     ",Edge.AdjacNode[ii])
print("\n\n")
print(len(Edge.AdjacElem))


#Attribute = np.zeros((Node.NNode,2))
#Attribute[:,0] = Node.u[0::2]
#Attribute[:,1] = Node.u[1::2]
#WriteResult.WriteCustomizedAttribute(Fem, Node, Element, Attribute)

toc = ct.time() - tic
Fem.PrintCommand(' ',0)
Fem.PrintCommand("Elapsed CPU time: "+str(toc)+"[sec]",0)
Fem.PrintCommand("file log is save in ./log/" +Fem.title + ".dat",0)
Fem.LogClose()

