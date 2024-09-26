# ----------------------------------------------------------------
# Written by: CK, HS in promechanics.org at Columbia University
# ----------------------------------------------------------------
import numpy as np
import sys
import time as ct
import os.path
from DataStructure import *

def GenerateEdge(Model):
    ModelFacet   = Model.Facet
    ModelElement = Model.Element
    ModelEdge    = []
    NEdge        = 0

    #print( "Edge start" )
    #print( "=================" )
    #print( "Dimension: ",Model.Dim )
    #print( "# of Elem: ",Model.NElem )
    #print( "# of Node: ",Model.NNode )
    #print( "# of Facet: ",Model.NFacet )
    #print( "=================" )

    for Facet in ModelFacet:
        facetNode = Facet.AdjNode.copy()
        facetNode.append(facetNode[0])

        for ii in range(len(facetNode)-1):
            flag = -1
            for jj, Edge in enumerate(ModelEdge):
                EdgeNode = Edge.AdjNode
                if EdgeNode[0] == facetNode[ii+1] and EdgeNode[1] == facetNode[ii]:
                    flag = jj
                    break
                elif EdgeNode[0] == facetNode[ii] and EdgeNode[1] == facetNode[ii+1]:
                    flag = jj
                    break

            if flag == -1:
                NEdge += 1
                AdjNode   = [facetNode[ii], facetNode[ii+1]]
                AdjFacet  = [Facet, Model.FacetInvalid]
                NewEdge   = EdgeAttribute(AdjNode, AdjFacet)
                ModelEdge.append(NewEdge)
            else:
                ModelEdge[flag].AdjFacet[1] = Facet

    # Example 1
    # print Node1 Node2 on Edge
    #for Edge in ModelEdge:
        #print(Edge.AdjNode[0].Id, Edge.AdjNode[1].Id)
    #print("===============================")

    # Example 2
    # print Facet1 - Edge - Facet2
    #for Edge in ModelEdge:
        #print(Edge.AdjNode[0].Id, Edge.AdjNode[1].Id)
        #print("------")
        #for ii in Edge.AdjFacet[0].AdjNode:
            #print(ii.Id)
        #print("------")
        #for ii in Edge.AdjFacet[1].AdjNode:
            #print(ii.Id)
        #print("===============================")

    print("\n\tTotal number of Edge: ",NEdge)
    return ModelEdge


def GenerateFacet(Model):
    ModelFacet   = []
    ModelElement = Model.Element
    NFacet       = 0
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

            for ii in range(len(elem)-1):
                flag = -1
                for jj, Facet in enumerate(ModelFacet):
                    # FacetNode: FacetNode
                    FacetNode = Facet.AdjNode
                    if FacetNode[0] == elem[ii+1] and FacetNode[1] == elem[ii]:
                        flag = jj
                        break

                if flag == -1:
                    NFacet += 1
                    AdjNode = [elem[ii], elem[ii+1]]
                    AdjElem = [Element, Model.ElementInvalid]
                    NewFacet  = FacetAttribute(AdjNode, AdjElem)
                    ModelFacet.append(NewFacet)
                else:
                    ModelFacet[flag].AdjElem[1] = Element

        #print("AdjNode")
        #for tmp in ModelFacet:
            #print(tmp.AdjNode[0].Id, tmp.AdjNode[1].Id)
        #print("AdjElem")
        #for tmp in ModelFacet:
            #print(tmp.AdjElem[0].Id, tmp.AdjElem[1].Id)
        #input(" ")
        #print("NFacet",NFacet)
        #exit(1)


    elif Model.Dim == 3:
        #assert False, "Dimension 3 is under construction"

        for Element in ModelElement:
            elem = Element.Connectivity.copy()
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

                for jj, Facet in enumerate(ModelFacet):
                    LocalFacet = list(reversed(Facet.AdjNode))
                    if LocalFacet in circular_shift:
                        flag = jj
                        break

                if flag == -1:
                    NFacet += 1
                    AdjNode = circular_shift[0]
                    AdjElem = [Element, Model.ElementInvalid]
                    NewFacet  = FacetAttribute(AdjNode, AdjElem)
                    ModelFacet.append(NewFacet)
                else:
                    ModelFacet[flag].AdjElem[1] = Element

    #print("AdjElem")
    #for ii in ModelFacet:
        ##print("----")
        #print(ii.AdjElem[0].Id, ii.AdjElem[1].Id, ii.AdjElem[0]==Model.ElementInvalid, ii.AdjElem[1]==Model.ElementInvalid)

    #print("AdjNode")
    #for ii in ModelFacet:
        ##print("----")
        #print(ii.AdjNode[0].Id, ii.AdjNode[1].Id, ii.AdjNode[2].Id, ii.AdjNode[3].Id)
    #exit(1)

    print( "\n\tTotal Number of Facet: ",NFacet )
    return ModelFacet


    #def GenerateFacet(Fem, Element):
        #Facet = FacetAttribute()
        #if Fem.Dimension == 2:
            #for kk in range(Element.NElem):
                #elem = Element.Connectivity[kk]
                #elem = elem.copy()
                #elem.append(elem[0])
                #for ii in range(len(elem)-1):
                    #flag = -1
                    #for jj, EN in enumerate(Facet.AdjNode):
                        #if EN[0] == elem[ii+1] and EN[1] == elem[ii]:
                            #flag = jj
                            #break

                    #if flag == -1:
                        #Facet.AdjNode.append([elem[ii], elem[ii+1]])
                        #Facet.AdjElem.append([Element.Id[kk],-1])
                    #else:
                        #Facet.AdjElem[flag][1] = Element.Id[kk]

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
                    #for jj, EN in enumerate(Facet.AdjNode):
                        #tmp = list(reversed(EN))
                        #if tmp in circular_shift:
                            #flag = jj
                            #break

                    #if flag == -1:
                        #Facet.AdjNode.append(circular_shift[0])
                        #Facet.AdjElem.append([Element.Id[kk],-1])
                    #else:
                        #Facet.AdjElem[flag][1] = Element.Id[kk]
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
            #for jj, EN in enumerate(Edge.AdjNode):
                #if EN[0] == elem[ii+1] and EN[1] == elem[ii]:
                    #flag = jj
                    #break

            #if flag == -1:
                #Edge.AdjNode.append([elem[ii], elem[ii+1]])
                #Edge.AdjElem.append([Element.Id[kk],-1])
            #else:
                #Edge.AdjElem[flag][1] = Element.Id[kk]
    #return Edge
