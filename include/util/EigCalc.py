
# Import necessary packages and functions
import numpy as np
from .coordinate_transforms import *

def so2_metric(truth, inp):
    err = np.sqrt((np.sin(truth)-np.sin(inp))**2 + (np.cos(truth)-np.cos(inp))**2 )
    return err

def LogSo3(matrix, theta):
    if abs(theta) < 1e-8:
        tmp = 0.
    else:
        tmp = theta / (2.0 * np.sin(theta))*(matrix - matrix.T)
    return tmp

def so3_metric(vc_truth, vc_inp):

    # Riemannian distance
    # https://lcvmwww.epfl.ch/publications/data/articles/63/simaxpaper.pdf
    # http://lucaballan.altervista.org/pdfs/IK.pdf
    # https://www.cs.cmu.edu/~cga/dynopt/readings/Rmetric.pdf
    #########################################
    #R = vc_truth @ vc_inp.T
    #traceR = np.trace(R)
    #if np.trace(R) <= -1.0:
        #traceR = -1.
    #theta = np.arccos((traceR-1.0)*0.5)
    #tmp = LogSo3(R, theta)
    #err = np.linalg.norm(tmp)

    #print("\nvc_truth\n",vc_truth)
    #print("\nvc_inp\n",vc_inp)
    #print("\nR\n",R)
    #print("\nnp.trace(R)\n",np.trace(R))
    #print("\ntheta\n",theta)
    #print("\nresult\n",tmp)
    #print("\nerr\n",err)
    #input("")
    #########################################

    # Euclidean distance
    err = np.linalg.norm(vc_truth - vc_inp)
    #print("\nresult\n",np.linalg.norm(tmp))

    # FYI
    # np.linalg.norm(tmp) = np.sqrt(tmp[0,0]**2 + tmp[0,1]**2 + tmp[0,2]**2\
                                #+ tmp[1,0]**2  + tmp[1,1]**2  + tmp[1,2]**2\
                                #+ tmp[2,0]**2  + tmp[2,1]**2  + tmp[2,2]**2)
    return err

#global tol

def GetNNInput(eigenvalue, eigenvector, theta_prev=0., alpha_prev=0., beta_prev=0., gamma_prev=0. ):

    ####np.set_printoptions(formatter={'float_kind': lambda x: "{0:0.3f}".format(x)})

    abg_measure   = []
    theta_measure = []

    alpha_list    = []
    beta_list     = []
    gamma_list    = []

    theta_list    = []

    combination_vl = []
    combination_vc = []
    vc_prev        = rotation_matrix(alpha_prev, beta_prev, gamma_prev)

    combination_vl.append([0,1,2])
    combination_vl.append([0,2,1])
    combination_vl.append([1,0,2])
    combination_vl.append([1,2,0])
    combination_vl.append([2,0,1])
    combination_vl.append([2,1,0])

    combination_vc.append([1,1,1])
    combination_vc.append([1,-1,-1])
    combination_vc.append([-1,1,-1])
    combination_vc.append([-1,-1,1])

    vc      = np.copy(eigenvector)
    vl      = np.copy(eigenvalue)
    vc_tmp  = np.zeros_like(eigenvector)
    vc_tmp1 = np.zeros_like(eigenvector)
    vl_tmp  = np.zeros_like(eigenvalue)

    ####################### theta #####################################
    #print("-"*45)
    for jj in combination_vl:
        vl_tmp = np.copy(vl[jj])
        p, rho, theta = convert_123_to_prt(vl_tmp[0], vl_tmp[1], vl_tmp[2])
        tmp1, tmp2, tmp3 = convert_prt_to_123(p, rho, theta)
        #print(p, rho, theta/np.pi*180.)
        #print(vl_tmp)
        #print(tmp1, tmp2, tmp3)
        #print("="*40)
        err = so2_metric(theta_prev, theta)
        theta_measure.append(err)
        theta_list.append(theta)
    #exit(1)
    #print("here")
    #print(theta_measure)
    #print(np.array(theta_list)/np.pi*180.0)
    #print(np.where(theta_measure< min(theta_measure)+0.1))
    ##########################################################################

    index1 = np.where(theta_measure < min(theta_measure)+0.1)
    for ii in index1[0]:
        vl_tmp = np.copy(vl[combination_vl[ii]]) # eigen values
        vc = eigenvector.T[combination_vl[ii]].T # eigen vectors
        vc[:,2] = np.cross(vc[:,0],vc[:,1])
        #print("--"*20)
        #print("\neigenvector\n",eigenvector)
        #print("\nvc\n",vc)
        #print("--"*20)

        x1 = np.copy(vc[:,0])
        x2 = np.copy(vc[:,1])
        x3 = np.copy(vc[:,2])
        if np.linalg.norm(np.cross(x1,x2)- x3) > 1e-4:
            print("\nEigenVector\n",vc)
            assert False, "Check the right hand rule"

        for combi in combination_vc:
            #print("@@@@"*10)
            #print(combi)
            for jj in range(3):
                vc_tmp[:,jj] = vc[:,jj]*combi[jj]
            #print("\nvc_tmp\n",vc_tmp)

            alp, beta, gamma = rotation_angles(vc_tmp)
            if np.abs(beta - np.pi/2)<1e-5:
                alp   = alpha_prev
                gamma -= alp
            elif np.abs(beta + np.pi/2)<1e-5:
                alp = alpha_prev
                gamma += alp
            err = so3_metric(vc_prev, vc_tmp)
            #print("\nerror\n",err)

            abg_measure.append( err )
            alpha_list.append(alp)
            beta_list.append(beta)
            gamma_list.append(gamma)
        #input("")

    index2 = np.argmin(abg_measure)
    #print("abg_measure\n",abg_measure)
    #print(index2)
    #print("="*50)

    p     = p
    rho   = rho
    theta = theta_list[index1[0][index2//4]]
    alpha = alpha_list[index2]
    beta  = beta_list[index2]
    gamma = gamma_list[index2]
    #eigenvalue  =
    #eigenvector = rotation_matrix(alpha,beta,gamma)

    return p, rho, theta, alpha, beta ,gamma

def SortEig2D(eigenvalue, eigenvector):
    n1 = np.copy(eigenvector[:,0])
    n2 = np.copy(eigenvector[:,1])

    error_list = []
    error_list.append(np.linalg.norm( n1 - np.array([1.,0.])))
    error_list.append(np.linalg.norm( n1 - np.array([-1.,0.])))
    error_list.append(np.linalg.norm( n2 - np.array([1.,0.])))
    error_list.append(np.linalg.norm( n2 - np.array([-1.,0.])))

    D = np.zeros((2))
    V = np.zeros((2,2))
    ind1 = np.argmin(np.array(error_list))
    if ind1 == 0 or ind1 == 1:
        V[:,0] = n1
        V[:,1] = n2
        D[0]   = eigenvalue[0]
        D[1]   = eigenvalue[1]
    else:
        V[:,0] = n2
        V[:,1] = n1
        D[0]   = eigenvalue[1]
        D[1]   = eigenvalue[0]
    return D, V




def SortEig(eigenvalue, eigenvector, alpha_prev=0., beta_prev=0., gamma_prev=0. ):

    n1 = eigenvector[:,0]
    n2 = eigenvector[:,1]
    n3 = np.cross(n1,n2)
    #n3 = eigenvector[:,2]

    vc_prev        = rotation_matrix(alpha_prev, beta_prev, gamma_prev)

    ################ Check ################
    #if np.linalg.norm(np.cross(n1,n2) + n3) < 1e-4:
        #n3 = np.cross(n1,n2)

    if np.linalg.norm(np.cross(n1,n2)- n3) > 1e-4:

        #print("\nn1",n1)
        #print("\nn2",n2)
        #print("\nn3",n3)
        #print(n1, "\t",n2, "\t",n3)
        print(eigenvector)
        assert False, "Check the right hand rule"
    ################ Check ###############

    combination_vl   = []
    combination_vc   = []
    ErrorMeasureList = []

    combination_vl.append([0,1,2])
    combination_vl.append([0,2,1])
    combination_vl.append([1,0,2])
    combination_vl.append([1,2,0])
    combination_vl.append([2,0,1])
    combination_vl.append([2,1,0])

    combination_vc.append([1,1,1])
    combination_vc.append([1,-1,-1])
    combination_vc.append([-1,1,-1])
    combination_vc.append([-1,-1,1])


    vc      = np.zeros((3,3))
    vc[:,0] = n1
    vc[:,1] = n2
    vc[:,2] = n3
    vl      = np.copy(eigenvalue)
    vc_tmp  = np.zeros_like(eigenvector)
    vc_tmp1 = np.zeros_like(eigenvector)
    vl_tmp  = np.zeros_like(eigenvalue)

    ####################### theta #####################################
    #print("-"*45)
    ind_tmp = 0
    for jj in combination_vl:
        for combi in combination_vc:
            #ind_tmp +=1
            vc_tmp = np.copy(vc.T[jj].T)
            vc_tmp[:,2] = np.cross(vc_tmp[:,0],vc_tmp[:,1])
            #print(vc_tmp)
            #print(np.cross(vc_tmp[:,0],vc_tmp[:,1]))

            vc_tmp[:,0] *= combi[0]
            vc_tmp[:,1] *= combi[1]
            vc_tmp[:,2] *= combi[2]
            #print(ind_tmp)
            #print(vc_tmp)
            #ind_tmp +=1
            #print("=================")
            err = so3_metric(vc_prev, vc_tmp)
            ErrorMeasureList.append(err)
        #print("-"*40)

    #print( ErrorMeasureList )

    ind1 = np.argmin(np.array(ErrorMeasureList))
    #for ii in ErrorMeasureList:
        #print(ii)
    #print(ErrorMeasureList)
    #print(ind1)
    #exit(1)
    ind2 = divmod(ind1 ,4)[0]
    ind3 = ind1 % 4
    D = np.copy(vl[combination_vl[ind2]])
    V = np.copy(vc.T[combination_vl[ind2]].T)
    V[:,2] = np.cross(V[:,0], V[:,1])
    combi   = combination_vc[ind3]
    V[:,0] *= combi[0]
    V[:,1] *= combi[1]
    V[:,2] *= combi[2]
    return D, V

def GetEulerAngle(eigenvalue, eigenvector, alpha_prev=0., beta_prev=0., gamma_prev=0. ):
    n1 = eigenvector[:,0]
    n2 = eigenvector[:,1]
    #n3 = eigenvector[:,2]
    n3 = np.cross(n1,n2)
    vc_standard = rotation_matrix(alpha_prev, beta_prev, gamma_prev)

    value = []
    vc    = np.zeros_like(eigenvector)
    vl    = np.zeros_like(eigenvalue)
    ii = 0
    alpha_list  = []
    beta_list   = []
    gamma_list  = []
    combination = []
    # Case1 #################################
    combination.append([0,1,2])
    vc[:,0] = n1
    vc[:,1] = n2
    vc[:,2] = n3
    x1 = np.copy(vc[:,0])
    x2 = np.copy(vc[:,1])
    x3 = np.copy(vc[:,2])
    if np.linalg.norm(np.cross(x1,x2)- x3) > 1e-4:
        print(x1,x2,x3)
        assert False, "Check the right hand rule"

    # print("="*40)
    # 1
    vc[:,0] = x1
    vc[:,1] = x2
    vc[:,2] = x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 2
    vc[:,0] = x1
    vc[:,1] = -x2
    vc[:,2] = -x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 3
    vc[:,0] = -x1
    vc[:,1] = x2
    vc[:,2] = -x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 4
    vc[:,0] = -x1
    vc[:,1] = -x2
    vc[:,2] = x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # Case2 #################################
    #vl[0]   = eigenvalue[0]
    #vl[1]   = eigenvalue[2]
    #vl[2]   = eigenvalue[1]
    combination.append([0,2,1])
    vc[:,0] = -n1
    vc[:,1] = n3
    vc[:,2] = n2
    x1 = np.copy(vc[:,0])
    x2 = np.copy(vc[:,1])
    x3 = np.copy(vc[:,2])
    if np.linalg.norm(np.cross(x1,x2)- x3) > 1e-4:
        print(x1,x2,x3)
        print(np.linalg.norm(np.cross(x1,x2)- x3))
        assert False, "Check the right hand rule"

    # print("="*40)
    # 1
    vc[:,0] = x1
    vc[:,1] = x2
    vc[:,2] = x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 2
    vc[:,0] = x1
    vc[:,1] = -x2
    vc[:,2] = -x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 3
    vc[:,0] = -x1
    vc[:,1] = x2
    vc[:,2] = -x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 4
    vc[:,0] = -x1
    vc[:,1] = -x2
    vc[:,2] = x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1


    # Case3 #################################
    #vl[0]   = eigenvalue[1]
    #vl[1]   = eigenvalue[0]
    #vl[2]   = eigenvalue[2]
    combination.append([1,0,2])
    vc[:,0] = n2
    vc[:,1] = n1
    vc[:,2] = -n3
    x1 = np.copy(vc[:,0])
    x2 = np.copy(vc[:,1])
    x3 = np.copy(vc[:,2])
    if np.linalg.norm(np.cross(x1,x2)- x3) > 1e-4:
        print(x1,x2,x3)
        assert False, "Check the right hand rule"


    # print("="*40)
    # 1
    vc[:,0] = x1
    vc[:,1] = x2
    vc[:,2] = x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 2
    vc[:,0] = x1
    vc[:,1] = -x2
    vc[:,2] = -x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 3
    vc[:,0] = -x1
    vc[:,1] = x2
    vc[:,2] = -x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 4
    vc[:,0] = -x1
    vc[:,1] = -x2
    vc[:,2] = x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1


    # Case4 #################################
    #vl[0]   = eigenvalue[1]
    #vl[1]   = eigenvalue[0]
    #vl[2]   = eigenvalue[2]
    combination.append([1,2,0])
    vc[:,0] = n2
    vc[:,1] = n3
    vc[:,2] = n1
    x1 = np.copy(vc[:,0])
    x2 = np.copy(vc[:,1])
    x3 = np.copy(vc[:,2])
    if np.linalg.norm(np.cross(x1,x2)- x3) > 1e-4:
        print(x1,x2,x3)
        assert False, "Check the right hand rule"


    # print("="*40)
    # 1
    vc[:,0] = x1
    vc[:,1] = x2
    vc[:,2] = x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 2
    vc[:,0] = x1
    vc[:,1] = -x2
    vc[:,2] = -x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 3
    vc[:,0] = -x1
    vc[:,1] = x2
    vc[:,2] = -x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 4
    vc[:,0] = -x1
    vc[:,1] = -x2
    vc[:,2] = x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1



    # Case5 #################################
    #vl[0]   = eigenvalue[2]
    #vl[1]   = eigenvalue[0]
    #vl[2]   = eigenvalue[1]
    combination.append([2,0,1])
    vc[:,0] = n3
    vc[:,1] = n1
    vc[:,2] = n2
    x1 = np.copy(vc[:,0])
    x2 = np.copy(vc[:,1])
    x3 = np.copy(vc[:,2])
    if np.linalg.norm(np.cross(x1,x2)- x3) > 1e-4:
        print(x1,x2,x3)
        assert False, "Check the right hand rule"

    # print("="*40)
    # 1
    vc[:,0] = x1
    vc[:,1] = x2
    vc[:,2] = x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 2
    vc[:,0] = x1
    vc[:,1] = -x2
    vc[:,2] = -x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 3
    vc[:,0] = -x1
    vc[:,1] = x2
    vc[:,2] = -x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 4
    vc[:,0] = -x1
    vc[:,1] = -x2
    vc[:,2] = x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1




    # Case6 #################################
    #vl[0]   = eigenvalue[2]
    #vl[1]   = eigenvalue[1]
    #vl[2]   = eigenvalue[0]
    combination.append([2,1,0])
    vc[:,0] = n3
    vc[:,1] = n2
    vc[:,2] = -n1
    x1 = np.copy(vc[:,0])
    x2 = np.copy(vc[:,1])
    x3 = np.copy(vc[:,2])
    if np.linalg.norm(np.cross(x1,x2)- x3) > 1e-4:
        print(x1,x2,x3)
        assert False, "Check the right hand rule"

    # print("="*40)
    # 1
    vc[:,0] = x1
    vc[:,1] = x2
    vc[:,2] = x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 2
    vc[:,0] = x1
    vc[:,1] = -x2
    vc[:,2] = -x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 3
    vc[:,0] = -x1
    vc[:,1] = x2
    vc[:,2] = -x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)
    ii += 1

    # 4
    vc[:,0] = -x1
    vc[:,1] = -x2
    vc[:,2] = x3
    alp, beta, gamma = rotation_angles(vc)
    if np.abs(beta - np.pi/2)<1e-5:
        alp   = alpha_prev
        gamma -= alp
    elif np.abs(beta + np.pi/2)<1e-5:
        alp = alpha_prev
        gamma += alp
    value.append(np.linalg.norm( vc-vc_standard ))
    alpha_list.append(alp)
    beta_list.append(beta)
    gamma_list.append(gamma)
    # print("-"*40,"\n",ii,"\nEigenVector\n",vc)
    # print("alpha,\tbeta,\tgamma\n",alp/np.pi*180.,"\t",beta/np.pi*180.,"\t",gamma/np.pi*180.)

    ind1 = np.argmin(np.array(value))
    ind2 = divmod(ind1 ,4)[0]
    ################### print ###################
    #print("\n\n","="*40)
    #print("norm values\n",value)
    #print(ind1)
    #print("alpha,\tbeta,\tgamma (degree)\n",alpha_list[ind1]*180.0/3.141592,"\t", beta_list[ind1]*180.0/3.141592,"\t", gamma_list[ind1]*180.0/3.141592)
    #print( rotation_matrix(alpha_list[ind1], beta_list[ind1], gamma_list[ind1]))
    ################### print ###################
    vl[0] = eigenvalue[combination[ind2][0]]
    vl[1] = eigenvalue[combination[ind2][1]]
    vl[2] = eigenvalue[combination[ind2][2]]

    return vl, alpha_list[ind1], beta_list[ind1], gamma_list[ind1]

# http://eecs.qmul.ac.uk/~gslabaugh/publications/euler.pdf
#Rx = [1 0 0
#      0 cos(alpha) -sin(alpha)
#      0 sin(alpha) cos(alpha)]

#Ry = [cos(beta)  0 sin(beta)
#       0         1    0
#      -sin(beta) 0 cos(beta)]

#Rz = [cos(gamma) -sin(gamma) 0
#      sin(gamma) cos(gamma) 0
#      0 0 1]
#The rotation matrix is Rz@Ry@Rx
#
#   R = [ cos(beta)*cos(gamma)      sin(alpha)*sin(beta)*cos(gamma)-cos(alpha)*sin(gamma)      cos(alpha)*sin(beta)*cos(gamma)+sin(alpha)*sin(gamma)
#         cos(beta)*sin(gamma)      sin(alpha)*sin(beta)*sin(gamma)+cos(alpha)*cos(gamma)      cos(alpha)*sin(beta)*sin(gamma)-sin(alpha)*cos(gamma)
#         -sin(beta)                sin(alpha)*cos(beta)                                       cos(alpha)*cos(beta)

#The rotation matrix is Rx@Ry@Rz
#
#   R = [ cos(beta)*cos(gamma)                                      -cos(beta)*sin(gamma)                                        sin(beta)
#         cos(alpha)*sin(gamma)+sin(alpha)*sin(beta)*cos(gamma)      cos(alpha)*cos(gamma)-sin(alpha)*sin(beta)*sin(gamma)      -sin(alpha)*cos(beta)
#         sin(alpha)*sin(gamma)-cos(alpha)*sin(beta)*cos(gamma)      sin(alpha)*cos(gamma)+cos(alpha)*sin(beta)*sin(gamma)      cos(alpha)*cos(beta)
#

def ReconMat(matrix):
    # This ensures that
    # alpha is between -pi and pi
    # beta  is between  0. and pi/2.
    # gamma is between  0. and pi
    sign1 = 1
    sign2 = 1
    if matrix[0,2] < 0.:
        sign2 = -1

    if matrix[0,1] > 0.:
        sign1 = -1

    matrix[:,2] = matrix[:,2] * sign2
    matrix[:,1] = matrix[:,1] * sign1
    matrix[:,0] = np.cross( matrix[:,1], matrix[:,2] )
    return matrix

def ReconMatGamma(A,R):
    vec  = np.array([0,0,1])
    tol1  = 1e-2
    tmp  = np.zeros((3,3))
    tmp1 = np.zeros((3,3))

    if np.linalg.norm(R[:,0] - vec) < tol1 or np.linalg.norm(R[:,0] + vec) < tol1:
        tmp[:,2] = R[:,0]
        tmp[:,0] = R[:,1]
        tmp[:,1] = R[:,2]
        tmp1[2,2] = A[0,0]
        tmp1[0,0] = A[1,1]
        tmp1[1,1] = A[2,2]

    elif np.linalg.norm(R[:,1] - vec) < tol1 or np.linalg.norm(R[:,1] + vec) < tol1:
        tmp[:,2] = R[:,1]
        tmp[:,0] = R[:,2]
        tmp[:,1] = R[:,0]
        tmp1[2,2] = A[1,1]
        tmp1[0,0] = A[2,2]
        tmp1[1,1] = A[0,0]

    elif np.linalg.norm(R[:,2] - vec) < tol1 or np.linalg.norm(R[:,2] + vec) < tol1:
        tmp[:,2] = R[:,2]
        tmp[:,0] = R[:,0]
        tmp[:,1] = R[:,1]
        tmp1[2,2] = A[2,2]
        tmp1[0,0] = A[0,0]
        tmp1[1,1] = A[1,1]
    else:
        print(R)
        print(np.linalg.norm(R[:,0] - vec), "\t",np.linalg.norm(R[:,0] + vec))
        print(np.linalg.norm(R[:,1] - vec), "\t",np.linalg.norm(R[:,1] + vec))
        print(np.linalg.norm(R[:,2] - vec), "\t",np.linalg.norm(R[:,2] + vec))
        assert False, "Problem on ReconMatGamma, 'Stress state is out of the plane' in EigCalc.py"

    R = tmp
    A = tmp1

    alp1, beta1, gamma1 = rotation_angles(tmp)
    if(np.abs(alp1)<tol1 and np.abs(beta1)<tol1 ):
        R = tmp
        return A, R

    tmp[:,0] = -R[:,0]
    tmp[:,1] = -R[:,1]
    tmp[:,2] = R[:,2]
    alp1, beta1, gamma1 = rotation_angles(tmp)
    if(np.abs(alp1)<tol1 and np.abs(beta1)<tol1 ):
        R = tmp
        return A, R

    tmp[:,0] = R[:,0]
    tmp[:,1] = -R[:,1]
    tmp[:,2] = -R[:,2]
    alp1, beta1, gamma1 = rotation_angles(tmp)
    if(np.abs(alp1)<tol1 and np.abs(beta1)<tol1 ):
        R = tmp
        return A, R

    tmp[:,0] = -R[:,0]
    tmp[:,1] = R[:,1]
    tmp[:,2] = -R[:,2]
    alp1, beta1, gamma1 = rotation_angles(tmp)
    if(np.abs(alp1)<tol1 and np.abs(beta1)<tol1 ):
        R = tmp
        return A, R


    print(alp1, beta1, gamma1)
    print("\n\n",A,"\n\n",R)
    assert False, "Problem on ReconMatGamma, 'Stress state is out of the plane' in EigCalc.py"


def rotation_matrix(theta1, theta2, theta3):
    c1 = np.cos(theta1)
    s1 = np.sin(theta1)
    c2 = np.cos(theta2)
    s2 = np.sin(theta2)
    c3 = np.cos(theta3)
    s3 = np.sin(theta3)

    #V = np.zeros((3,3))

    #V[0,0] = c2*c3
    #V[0,1] =-c2*s3
    #V[0,2] = s2

    #V[1,0] = c1*s3+s1*s2*c3
    #V[1,1] = c1*c3-s1*s2*s3
    #V[1,2] =-c2*s1

    #V[2,0] = s1*s3-c1*s2*c3
    #V[2,1] = s1*c3+c1*s2*s3
    #V[2,2] = c1*c2
    #print("\nV\n",V)

    Ra = np.array([[1., 0., 0.],
                  [ 0., c1,-s1],
                  [ 0., s1, c1]])

    Rb = np.array([[ c2, 0., s2],
                   [ 0., 1., 0.],
                   [-s2, 0., c2]])

    Rc = np.array([[c3,-s3, 0.],
                   [s3, c3, 0.],
                   [0., 0., 1.]])

    matrix = Ra @ Rb @ Rc

    return matrix



def rotation_angles(matrix):
#The rotation matrix is Rx@Ry@Rz
#
#   R = [ cos(beta)*cos(gamma)                                      -cos(beta)*sin(gamma)                                        sin(beta)
#         cos(alpha)*sin(gamma)+sin(alpha)*sin(beta)*cos(gamma)      cos(alpha)*cos(gamma)-sin(alpha)*sin(beta)*sin(gamma)      -sin(alpha)*cos(beta)
#         sin(alpha)*sin(gamma)-cos(alpha)*sin(beta)*cos(gamma)      sin(alpha)*cos(gamma)+cos(alpha)*sin(beta)*sin(gamma)      cos(alpha)*cos(beta)
#
    # Rx@Ry@Rz
    # range of arctan2 is between -pi    and pi
    # range or arcsin  is between -pi/2. and pi/2.
    tol =1e-7

    r11, r12, r13 = matrix[0]
    r21, r22, r23 = matrix[1]
    r31, r32, r33 = matrix[2]

    beta = np.arcsin(r13)
    alp  = np.arctan2(-r23*np.cos(beta), r33*np.cos(beta))
    gam  = np.arctan2(-r12*np.cos(beta), r11*np.cos(beta))
    if (abs((beta - np.pi*0.5))<1e-5 or abs((beta + np.pi*0.5))<1e-5):
        alp = 0
        gam = np.arctan2(r21, r22)

    return alp, beta, gam




def CalcEig(A):
    tol= 1e-8
    R = np.eye(3,3)
    ind = 0
    A_org = np.copy(A)

    if np.max(np.abs(A)) < tol:
        return A, R

    while( np.linalg.norm( (A/np.max(np.abs(A)) - np.diag(np.diag(A/np.max(np.abs(A)))) )/np.sqrt(2.) ) > tol):
        ind += 1
        if( np.abs(A[1,2])/np.max(np.abs(A_org))>tol):
            #print("====Z\n",A)
            A, R = RotateX(A,R)
            #print("\n",A)

        if( np.abs(A[0,2])/np.max(np.abs(A_org))>tol):
            #print("====Y\n",A)
            A, R = RotateY(A,R)
            #print("\n",A)
            #input("")

        if( np.abs(A[0,1])/np.max(np.abs(A_org))>tol):
            #print("====Z\n",A)
            A, R = RotateZ(A,R)
            #print("\n",A)

        if( ind >2000):
            print("crt", np.linalg.norm( (A/np.max(np.abs(A)) - np.diag(np.diag(A/np.max(np.abs(A)))))/np.sqrt(2.) ) )
            print(ind)
            print("\nOriginal matrix\n",A_org)
            print("\nEigenValues\n",A)
            print("\nEigenVectors\n",R)
            vc,vl = np.linalg.eig(A_org)
            print("numpy calculation\n",vc,"\n",vl)
            input("eigenvalue calculation takes more ",ind," iterations")
    #print(ind)
    R = ReconMat(R)
    return A, R

def RotateX(A,R):
    alpha = np.arctan2(  2*A[1,2], A[2,2]-A[1,1] ) * 0.5

    Ra = np.array([[1., 0.,         0.],
                   [0., np.cos(alpha), -np.sin(alpha)],
                   [0., np.sin(alpha), np.cos(alpha)]])

    A  = Ra @ A @ np.transpose(Ra)
    R = R @ np.transpose(Ra)
    return A, R

def RotateY(A,R):
    beta = np.arctan2(  2*A[0,2], A[0,0]-A[2,2] ) * 0.5

    Rb = np.array([[np.cos(beta),  0, np.sin(beta)],
                   [0,         1,    0],
                   [-np.sin(beta), 0, np.cos(beta)]])

    A  = Rb @ A @ np.transpose(Rb)
    R = R @ np.transpose(Rb)
    return A, R

def RotateZ(A,R):
    gamma = np.arctan2(  2*A[0,1], A[1,1]-A[0,0] ) * 0.5

    Rc = np.array([[np.cos(gamma), -np.sin(gamma), 0],
                   [np.sin(gamma), np.cos(gamma), 0],
                   [0, 0, 1]])

    A  = Rc @ A @ np.transpose(Rc)
    R = R @ np.transpose(Rc)
    return A, R


def CalcEig_HS(sigma):
    tol =1e-6
    #sigma = np.array([[sig11, sig12, sig31],
                    #[sig12, sig22, sig23],
                    #[sig31, sig23, sig33]])

    D, V = np.linalg.eig(sigma)
    mask = np.abs(V/np.abs(np.max(V))) < tol
    V[mask] = 0.

    # Hyoung Suk Suh
    n1_tmp = V[:,0]
    n2_tmp = V[:,1]
    n3_tmp = V[:,2]

    evaluated_distances = np.ones(3)
    for i in range(np.shape(evaluated_distances)[0]):
        if i == 0:
            evaluated_distances[i] = np.abs(np.dot(n1_tmp, np.array([0.,0.,1.])))
        elif i == 1:
            evaluated_distances[i] = np.abs(np.dot(n2_tmp, np.array([0.,0.,1.])))
        elif i == 2:
            evaluated_distances[i] = np.abs(np.dot(n3_tmp, np.array([0.,0.,1.])))

    max_idx = np.argmax(evaluated_distances)

    if max_idx == 0: # SWAP 1 and 3
        sig_principal_mag = np.zeros_like(D)
        sig_principal_mag[0] = D[2]
        sig_principal_mag[1] = D[1]
        sig_principal_mag[2] = D[0]

        n1 = np.zeros_like(n1_tmp)
        n2 = np.zeros_like(n1_tmp)
        n3 = np.zeros_like(n1_tmp)

        n1[0] = n3_tmp[0]
        n1[1] = n3_tmp[1]
        n1[2] = n3_tmp[2]

        n3[0] = n1_tmp[0]
        n3[1] = n1_tmp[1]
        n3[2] = n1_tmp[2]

        n2 = np.cross(n3, n1)

        sig_principal_vec = np.zeros_like(V)
        sig_principal_vec[0,0] = n1[0]
        sig_principal_vec[1,0] = n1[1]
        sig_principal_vec[2,0] = n1[2]
        sig_principal_vec[0,1] = n2[0]
        sig_principal_vec[1,1] = n2[1]
        sig_principal_vec[2,1] = n2[2]
        sig_principal_vec[0,2] = n3[0]
        sig_principal_vec[1,2] = n3[1]
        sig_principal_vec[2,2] = n3[2]

    elif max_idx == 1: # SWAP 2 and 3
        sig_principal_mag = np.zeros_like(D)
        sig_principal_mag[0] = D[0]
        sig_principal_mag[1] = D[2]
        sig_principal_mag[2] = D[1]

        n1 = np.zeros_like(n1_tmp)
        n2 = np.zeros_like(n1_tmp)
        n3 = np.zeros_like(n1_tmp)

        n2[0] = n3_tmp[0]
        n2[1] = n3_tmp[1]
        n2[2] = n3_tmp[2]

        n3[0] = n2_tmp[0]
        n3[1] = n2_tmp[1]
        n3[2] = n2_tmp[2]

        n1 = np.cross(n2, n3)

        sig_principal_vec = np.zeros_like(V)
        sig_principal_vec[0,0] = n1[0]
        sig_principal_vec[1,0] = n1[1]
        sig_principal_vec[2,0] = n1[2]
        sig_principal_vec[0,1] = n2[0]
        sig_principal_vec[1,1] = n2[1]
        sig_principal_vec[2,1] = n2[2]
        sig_principal_vec[0,2] = n3[0]
        sig_principal_vec[1,2] = n3[1]
        sig_principal_vec[2,2] = n3[2]

    else:
        sig_principal_mag = D
        sig_principal_vec = V


    D = sig_principal_mag
    V = sig_principal_vec

    ###
    if V[2,2] < 0:
        V = -V

    n1_tmp = V[:,0]
    n2_tmp = V[:,1]
    n3_tmp = V[:,2]

    evaluated_distances = np.ones(3)
    for i in range(np.shape(evaluated_distances)[0]):
        if i == 0:
            evaluated_distances[i] = np.abs(np.dot(n1_tmp, np.array([1.,0.,0.])))
        elif i == 1:
            evaluated_distances[i] = np.abs(np.dot(n2_tmp, np.array([1.,0.,0.])))
        elif i == 2:
            evaluated_distances[i] = np.abs(np.dot(n3_tmp, np.array([1.,0.,0.])))

    if max_idx == 1: # SWAP 1 & 2
        sig_principal_mag = np.zeros_like(D)
        sig_principal_mag[0] = D[1]
        sig_principal_mag[1] = D[0]
        sig_principal_mag[2] = D[2]

        n1 = np.zeros_like(n1_tmp)
        n2 = np.zeros_like(n1_tmp)
        n3 = np.zeros_like(n1_tmp)

        n2[0] = n1_tmp[0]
        n2[1] = n1_tmp[1]
        n2[2] = n1_tmp[2]

        n3[0] = n3_tmp[0]
        n3[1] = n3_tmp[1]
        n3[2] = n3_tmp[2]

        n1 = np.cross(n2, n3)

        sig_principal_vec = np.zeros_like(V)
        sig_principal_vec[0,0] = n1[0]
        sig_principal_vec[1,0] = n1[1]
        sig_principal_vec[2,0] = n1[2]
        sig_principal_vec[0,1] = n2[0]
        sig_principal_vec[1,1] = n2[1]
        sig_principal_vec[2,1] = n2[2]
        sig_principal_vec[0,2] = n3[0]
        sig_principal_vec[1,2] = n3[1]
        sig_principal_vec[2,2] = n3[2]

        D = sig_principal_mag
        V = sig_principal_vec
    ###

    if V[0,0] < 0:
        V[:,0] = -V[:,0]
        V[:,1] = -V[:,1]

    if V[0,1] < 0:
        V[0,1] = -V[0,1]
        V[1,0] = -V[1,0]

    return D, V
