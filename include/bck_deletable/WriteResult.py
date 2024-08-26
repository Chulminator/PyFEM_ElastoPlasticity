# ---------------------------------------------------------------- 
# Written by: CK, HS in promechanics.org at Columbia University
# ----------------------------------------------------------------       
import numpy as np
import sys
import time as ct
import os.path

# Things to upgrade
#   PlotSetup   : what variable will be plotted how often
#   OutputSetup : what variable will be written how often


def WriteDisp(Fem, Node, Element):

    fid1 = open(Fem.result+'/'+Fem.title+'.NODE', "w")
    fid1.write("%s\n" % (Node.NNode))
    DispX = Node.u[0::2]
    DispY = Node.u[1::2]
    for ii in range(Node.NNode):
        fid1.write("%s\t%s\t%s\t%s\t%s\n" %(Node.Id[ii], Node.Coord[ii][0], Node.Coord[ii][1], DispX[ii], DispY[ii]))
    fid1.close()

    fid2 = open(Fem.result+'/'+Fem.title+'.ELEM', "w")
    fid2.write("%s\n" % (Element.NElem))
    for ii in range(Element.NElem):
        #print(Element.Id[ii])
        #print(type(Element.Id[ii]))
        #print("here\n\n")
        #fid2.write("%s\t", %(Element.Id[ii]) )
        fid2.write("%s\t" % (Element.Id[ii]))
        for jj in Element.Connectivity[ii]:
            fid2.write("%s\t" % (jj))
            #fid2.write("%s\t", %( jj ))
        #fid2.write("-1\n")
        fid2.write("\n")
    fid2.close()

def WriteStress(Fem, Node, Element):
    fid1 = open(Fem.result+'/'+Fem.title+'_Stress.NODE', "w")
    fid1.write("%s\n" % (Node.NNode))
    DispX = Node.u[0::2]
    DispY = Node.u[1::2]
    for ii in range(Node.NNode):
        fid1.write("%s\t%s\t%s\t%s\t%s\t%s\n" %(Node.Id[ii], Node.Coord[ii][0], Node.Coord[ii][1], Node.stress[ii,0], Node.stress[ii,1], Node.stress[ii,2]))
    fid1.close()

    fid2 = open(Fem.result+'/'+Fem.title+'_Stress.ELEM', "w")
    fid2.write("%s\n" % (Element.NElem))
    for ii in range(Element.NElem):
        #print(Element.Id[ii])
        #print(type(Element.Id[ii]))
        #print("here\n\n")
        #fid2.write("%s\t", %(Element.Id[ii]) )
        fid2.write("%s\t" % (Element.Id[ii]))
        for jj in Element.Connectivity[ii]:
            fid2.write("%s\t" % (jj))
        #fid2.write("-1\n")
        fid2.write("\n")
    fid2.close()


def WriteCustomizedAttribute(Fem, Node, Element, att):
    fid1 = open(Fem.result+'/'+Fem.title+'_CustomizedAttribute.NODE', "w")
    fid1.write("%s\n" % (Node.NNode))
    DispX = Node.u[0::2]
    DispY = Node.u[1::2]
    if not(att.shape[0] == Node.NNode):
        print( "The attribute size should be [NNode, NAttribute]; tmp>1" )
        print( "Check att.shape[0]" )

    if (att.shape[1] >8 ):
        print( "Warning: Matlab might not be able to read # of Attributes > 8 " )

    for ii in range(Node.NNode):
        fid1.write("%s\t%s\t%s\t" %(Node.Id[ii], Node.Coord[ii][0], Node.Coord[ii][1]))
        for jj in range(att.shape[1]):
            fid1.write("%s\t" %(att[ii][jj]))
        fid1.write("\n")
    fid1.close()

    fid2 = open(Fem.result+'/'+Fem.title+'_CustomizedAttribute.ELEM', "w")
    fid2.write("%s\n" % (Element.NElem))
    for ii in range(Element.NElem):
        fid2.write("%s\t" % (Element.Id[ii]))
        for jj in Element.Connectivity[ii]:
            fid2.write("%s\t" % (jj))
        #fid2.write("-1\n")
        fid2.write("\n")
    fid2.close()
