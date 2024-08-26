 

import pygmsh
import meshio
import numpy as np

def WriteNodeAndElementFromSandiaMesh(directory, Node, Element, Dim=2):


    NNode = Node.shape[0]
    NElem = Element.shape[0]
    fid1 = open(directory+'_'+str(NNode)+'_'+str(NElem)+'.NODE', "w")
    fid1.write("%s\n" % (len(Node)))

    for ii, node in enumerate(Node):
        if Dim == 2:
            fid1.write("%s\t%s\t%s\n" %(ii+1, node[0], node[1]))
        elif Dim == 3:
            fid1.write("%s\t%s\t%s\t%s\n" %(ii+1, node[0], node[1], node[2]))
        else:
            assert false, "Check the 4th input 'Dim'"


    fid1.close()

    fid2 = open(directory+'_'+str(NNode)+'_'+str(NElem)+'.ELEM', "w")
    fid2.write("%s\n" % (len(Element)))
    for ii, Connectivity in enumerate(Element):
        fid2.write("%s\t" % (ii+1))
        for jj in Connectivity:
            fid2.write("%s\t" % (jj+1))
        fid2.write("\n")
    fid2.close()

    return

def WriteNodeAndElementFromMeshio(mesh, directory, MeshType):

    fid1 = open(directory+'.NODE', "w")
    fid1.write("%s\n" % (len(mesh.points)))
    Node = []
    for ii, node in enumerate(mesh.points):
        fid1.write("%s\t%s\t%s\n" %(ii+1, node[0], node[1]))
        tmp =[]
        tmp.append(ii+1)
        tmp.append(node[0])
        tmp.append(node[1])
        Node.append(tmp)
    fid1.close()

    fid2 = open(directory+'.ELEM', "w")
    Element =[]
    if MeshType == "q4":
        fid2.write("%s\n" % (len(mesh.cells_dict["quad"])))
        for ii, Connectivity in enumerate(mesh.cells_dict["quad"]):
            fid2.write("%s\t" % (ii+1))
            tmp = []
            tmp.append(ii+1)
            for jj in Connectivity:
                fid2.write("%s\t" % (jj+1))
                tmp.append(jj+1)
            fid2.write("\n")
            Element.append(tmp)
        Node.append(tmp)
    fid2.close()
    return Node, Element

def MakeMeshioFromNodeAndElement(Node, Element, Fem):
    print("x")
    exit(1)
    return mesh
    #Mesh
