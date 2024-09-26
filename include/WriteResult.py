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

def WritePVD(Program, Model):
    name = Program.result+'/'+Program.title+'.pvd'
    fid1 = open(name, "w")
    fid1.write("<?xml version=\"1.0\"?>\n")
    fid1.write("<VTKFile type=\"Collection\" version=\"0.1\">\n")
    fid1.write("  <Collection>\n")
    for ii in range(Model.totalstep):
        tmp = (ii+1)/Model.totalstep
        fid1.write("    <DataSet timestep=\""+str(tmp)+"\" part=\"0\" file=\""+Program.title+"_"+str(ii+1)+".vtu\" />\n")

    fid1.write("  </Collection>\n")
    fid1.write("</VTKFile>\n")
    fid1.close()

def WriteVTU(Program, Model, att=None):
    # https://github.com/nschloe/meshio/blob/main/src/meshio/_common.py
    # https://pypi.org/project/pygmsh/4.0.8/#description
    #https://pypi.org/project/meshio/
    #import pygmsh
    import meshio


    #if Fem.Dimension == 3:
        #points = Node.Coord # for the plot with displacement, this should be changed
    #elif Fem.Dimension == 2:
        #for coord in Node.Coord:
            #coord.append(0.0)
        #points = Node.Coord
        #print(Node.Coord)
        #exit(1)
    #points = Node.Coord
    Connectivity = []
    points       = []

    for element in Model.Element:
        connec = []
        for node in element.Connectivity:
            connec.append(Model.Node.index(node))   # need to get checked
        Connectivity.append(connec)

    for node in Model.Node:
        points.append(node.Coord + node.u)

    if Model.Element[0].ElmentType == 'hex8':
        cells = [("hexahedron", Connectivity)]

    elif Model.Element[0].ElmentType == 'tet4':
        cells = [("tetra", Connectivity)]

    elif Model.Element[0].ElmentType == 't3':
        cells = [("triangle", Connectivity)]

    elif Model.Element[0].ElmentType == 'q4':
        cells = [("quad", Connectivity)]

    elif Model.Element[0].ElmentType == 't6':
        print("Model.ElmentType: ",Model.Element[0].ElmentType)
        assert False, "This has to be checked"
        cells = [("triangle", Connectivity)]

    else:
        print("Model.ElmentType",Model.Element[0].ElmentType)
        assert False, "Element type:  is not ready"

    if Model.Dim == 3:
        mesh = meshio.Mesh(
            points,
            cells,
            # Optionally provide extra data on points, cells, etc.
            #point_data={"VonMises": att[:,-2]},
            point_data={"VonMises": Model.NodalsigmaVM,
                        "lambda": Model.NodalPlasticMultiplier,
                        "dx": Model.Nodalu[0::3],
                        "dy": Model.Nodalu[1::3],
                        "dz": Model.Nodalu[2::3],
                        "e11": Model.Nodalstrain_e[:,0]+Model.Nodalstrain_p[:,0],
                        "e22": Model.Nodalstrain_e[:,1]+Model.Nodalstrain_p[:,1],
                        "e33": Model.Nodalstrain_e[:,2]+Model.Nodalstrain_p[:,2],
                        "e12": Model.Nodalstrain_e[:,3]+Model.Nodalstrain_p[:,3],
                        "e23": Model.Nodalstrain_e[:,4]+Model.Nodalstrain_p[:,4],
                        "e31": Model.Nodalstrain_e[:,5]+Model.Nodalstrain_p[:,5],
                        "s11": Model.Nodalstress[:,0],
                        "s22": Model.Nodalstress[:,1],
                        "s33": Model.Nodalstress[:,2],
                        "s12": Model.Nodalstress[:,3],
                        "s23": Model.Nodalstress[:,4],
                        "s31": Model.Nodalstress[:,5]},

            # Each item in cell data must match the cells array
            #cell_data={"a": [[0.1, 0.2], [0.4]]},
        )
    elif Model.Dim == 2:
        #if Fem.twoD == "planestress":
        mesh = meshio.Mesh(
            points,
            cells,
            # Optionally provide extra data on points, cells, etc.
            #point_data={"VonMises": att[:,-2]},
            point_data={"VonMises": Model.NodalsigmaVM,
                        "lambda": Model.NodalPlasticMultiplier,
                        "dx": Model.Nodalu[0::2],
                        "dy": Model.Nodalu[1::2],
                        "e11": Model.Nodalstrain_e[:,0]+Model.Nodalstrain_p[:,0],
                        "e22": Model.Nodalstrain_e[:,1]+Model.Nodalstrain_p[:,1],
                        "e12": Model.Nodalstrain_e[:,2]+Model.Nodalstrain_p[:,2],
                        "s11": Model.Nodalstress[:,0],
                        "s22": Model.Nodalstress[:,1],
                        "s12": Model.Nodalstress[:,2]},

            # Each item in cell data must match the cells array
            #cell_data={"a": [[0.1, 0.2], [0.4]]},
        )
        #elif Fem.twoD == "planestrain":
            #mesh = meshio.Mesh(
                #points,
                #cells,
                ## Optionally provide extra data on points, cells, etc.
                ##point_data={"VonMises": att[:,-2]},
                #point_data={"VonMises": Node.sigmaVM,
                            #"lambda": Node.PlasticMultiplier,
                            #"dx": Node.u[0::2],
                            #"dy": Node.u[1::2],
                            #"e11": Node.strain_e[:,0]+Node.strain_p[:,0],
                            #"e22": Node.strain_e[:,1]+Node.strain_p[:,1],
                            #"e12": Node.strain_e[:,2]+Node.strain_p[:,2],
                            #"s11": Node.stress[:,0],
                            #"s22": Node.stress[:,1],
                            #"s12": Node.stress[:,2]},

                ## Each item in cell data must match the cells array
                ##cell_data={"a": [[0.1, 0.2], [0.4]]},
            #)
            #s
        #else:
            #assert False, "Check Fem.twoD"





    mesh.write(
        Program.result+'/'+Program.title+'_'+str(Model.step)+".vtu",
    )
    return



def WriteCustomizedAttribute(Program, Model, att): # edit
    fid1 = open(Program.result+'/'+Program.title+'_CustomizedAttribute_'+str(Model.step)+'.NODE', "w")
    fid1.write("%s\n" % (Model.NNode))


    if not(att.shape[0] == Model.NNode):
        print( "The attribute size should be [NNode, NAttribute]; tmp>1" )
        print( "Check att.shape[0]" )

    if (att.shape[1] >8 ):
        print( "WARNING: Matlab might not be able to read # of Attributes > 8 " )

    for ii, node in enumerate(Model.Node):
        fid1.write("%s\t" %(node.Id))
        for jj in node.Coord:
            fid1.write("%s\t" %(jj))
        for jj in att[ii]:
            fid1.write("%s\t" %(jj))
        fid1.write("\n")
    fid1.close()

    fid2 = open(Program.result+'/'+Program.title+'_CustomizedAttribute_'+str(Model.step)+'.ELEM', "w")
    fid2.write("%s\n" % (Model.NElem))
    for ii, elem in enumerate(Model.Element):
        fid2.write("%s\t" % (elem.Id))
        for node in elem.Connectivity:
            fid2.write("%s\t" % (node.Id))
        fid2.write("\n")
    fid2.close()
