# ----------------------------------------------------------------
# Written by: CK, HS in promechanics.org at Columbia University
# ----------------------------------------------------------------
import numpy as np
import sys
import time as ct
import os.path
from DataStructure import *

#class EdgeAttribute:
    #def __init__(self, Model):
        #self.AdjacFacet = []
        #self.AdjacNode  = []
        #self.NEdge      = 0

        #ModelFacet = Model.Facet
        #print( "Facet start" )
        #print( "=================" )
        #print( "Dimension: ",Model.Dim )
        #print( "# of Elem: ",Model.NElem )
        #print( "# of Node: ",Model.NNode )
        #print( "# of Facet: ",Model.Facet )
        #print( "=================" )
        #if Model.Dim == 2:
            #for Facet in ModelFacet:
                #exit(1)
        #elif Model.Dim == 3:
            #for LocalFacet in ModelFacet:
                #for Facet in LocalFacet.AdjacNode:
                    #print(Facet.Id)
                #exit(1)
                ##Facet = Element.Connectivity.copy()



        #exit(1)


class FacetAttribute:
    def __init__(self, Model):
        self.AdjacElem = []
        self.AdjacNode = []
        self.NFacet    = 0

        ModelElement = Model.Element
        #print( "Facet start" )
        #print( "=================" )
        #print( "Dimension: ",Model.Dim )
        #print( "# of Elem: ",Model.NElem )
        #print( "# of Node: ",Model.NNode )
        #print( "=================" )
        #exit(1)
        if Model.Dim == 2:
            for Element in ModelElement:
                elem = Element.Connectivity.copy()
                elem.append(elem[0])
                #for tmp in elem:
                    #print(tmp.Id)
                #exit(1)

                for ii in range(len(elem)-1):
                    flag = -1
                    for jj, EN in enumerate(self.AdjacNode):
                        # EN: FacetNode
                        if EN[0] == elem[ii+1] and EN[1] == elem[ii]:
                            flag = jj
                            break

                    if flag == -1:
                        self.NFacet += 1
                        self.AdjacNode.append([elem[ii], elem[ii+1]])
                        self.AdjacElem.append([Element, Model.ElementInvalid])

                    else:
                        self.AdjacElem[flag][1] = Element

                #print("AdjacNode")
                #for tmp in self.AdjacNode:
                    #print(tmp[0].Id, tmp[1].Id)
                #print("AdjacElem")
                #for tmp in self.AdjacElem:
                    #print(tmp[0].Id, tmp[1].Id)
                #input(" ")
            #print("NFacet",self.NFacet)
            #exit(1)


        elif Model.Dim == 3:
            #assert False, "Dimension 3 is under construction"

            for Element in ModelElement:
                elem = Element.Connectivity.copy()

                print( "Element id : ", Element.Id )
                for tmp in elem:
                    print(tmp.Id)
                print( "---------" )
                #exit(1)


                if Element.ElmentType == 'hex8':
                    facet_list = [[4,5,6,7],
                                  [5,1,2,6],
                                  [1,0,3,2],
                                  [0,4,7,3],
                                  [7,6,2,3],
                                  [0,1,5,4]]
                elif Element.ElmentType == 'tet4':
                    print("Warning: Tet4 is not fully checked!")
                    input("Press enter to proceed the analysis")
                    #assert False, "Element.ElmentType not ready"
                    facet_list = [[0,1,2],
                                  [0,1,3],
                                  [1,2,3],
                                  [2,0,3]]
                else:
                    assert False, "Element.ElmentType not ready"

                for LFacet in facet_list:
                    # LFacet: Local facet
                    # GFacet: Global facet

                    GFacet = []
                    for ii in LFacet:
                        GFacet.append( elem[ii] )

                    flag = -1
                    circular_shift = []
                    circular_shift.append(GFacet.copy())
                    for jj in range(len(GFacet)-1):
                        GFacet.append(GFacet.pop(0))
                        tmp = GFacet.copy()
                        circular_shift.append(tmp)

                    #print("=======")
                    #for jj in circular_shift:
                        #for tmp in jj:
                            #print( tmp.Id)
                        #print("------")
                    #exit(1)
                    #print("LFacet\t",LFacet)

                    #print("GFacet\t",np.array(GFacet)-1)

                    for jj, EN in enumerate(self.AdjacNode):
                        tmp = list(reversed(EN))
                        if tmp in circular_shift:
                            flag = jj
                            break

                    if flag == -1:
                        self.NFacet += 1
                        self.AdjacNode.append(circular_shift[0])
                        self.AdjacElem.append([Element,Model.ElementInvalid])
                    else:
                        self.AdjacElem[flag][1] = Element
                #print("")
            #print("AdjacElem")
            #for ii in self.AdjacElem:
                ##print("----")
                #print(ii[0].Id, ii[1].Id)

            #print("AdjacNode")
            #for ii in self.AdjacNode:
                ##print("----")
                #print(ii[0].Id, ii[1].Id, ii[2].Id, ii[3].Id)
            #exit(1)
        return


    #def GenerateFacet(Fem, Element):
        #Facet = FacetAttribute()
        #if Fem.Dimension == 2:
            #for kk in range(Element.NElem):
                #elem = Element.Connectivity[kk]
                #elem = elem.copy()
                #elem.append(elem[0])
                #for ii in range(len(elem)-1):
                    #flag = -1
                    #for jj, EN in enumerate(Facet.AdjacNode):
                        #if EN[0] == elem[ii+1] and EN[1] == elem[ii]:
                            #flag = jj
                            #break

                    #if flag == -1:
                        #Facet.AdjacNode.append([elem[ii], elem[ii+1]])
                        #Facet.AdjacElem.append([Element.Id[kk],-1])
                    #else:
                        #Facet.AdjacElem[flag][1] = Element.Id[kk]

        #elif Fem.Dimension==3:
            #for kk in range(Element.NElem):
                #elem = Element.Connectivity[kk]
                #elem = np.array(elem.copy())
                #if Fem.ElmentType == 'hex8':
                    #facet_list = [[4,5,6,7],
                                #[5,1,2,6],
                                #[1,0,3,2],
                                #[0,4,7,3],
                                #[7,6,2,3],
                                #[0,1,5,4]]
                #elif Fem.ElmentType == 'tet4':
                    #assert False, "Fem.ElmentType not ready"
                #else:
                    #assert False, "Fem.ElmentType not ready"
                #for LFacet in facet_list:
                    #GFacet = list(elem[LFacet])
                    #flag = -1
                    #circular_shift = []
                    #circular_shift.append(GFacet.copy())
                    #for jj in range(len(GFacet)-1):
                        #GFacet.append(GFacet.pop(0))
                        #tmp = GFacet.copy()
                        #circular_shift.append(tmp)

                    ##print("LFacet\t",LFacet)
                    #print("GFacet\t",np.array(GFacet)-1)
                    #for jj, EN in enumerate(Facet.AdjacNode):
                        #tmp = list(reversed(EN))
                        #if tmp in circular_shift:
                            #flag = jj
                            #break

                    #if flag == -1:
                        #Facet.AdjacNode.append(circular_shift[0])
                        #Facet.AdjacElem.append([Element.Id[kk],-1])
                    #else:
                        #Facet.AdjacElem[flag][1] = Element.Id[kk]
                #print("")

        #return Facet


#def GenerateEdge(Element):
    #Edge = EdgeAttribute()
    #for kk in range(Element.NElem):
        ##print()
        #elem = Element.Connectivity[kk]
        #elem = elem.copy()
        #elem.append(elem[0])
        #for ii in range(len(elem)-1):
            #flag = -1
            #for jj, EN in enumerate(Edge.AdjacNode):
                #if EN[0] == elem[ii+1] and EN[1] == elem[ii]:
                    #flag = jj
                    #break

            #if flag == -1:
                #Edge.AdjacNode.append([elem[ii], elem[ii+1]])
                #Edge.AdjacElem.append([Element.Id[kk],-1])
            #else:
                #Edge.AdjacElem[flag][1] = Element.Id[kk]
    #return Edge
