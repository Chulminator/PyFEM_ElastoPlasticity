# Generated with SMOP  0.41
from libsmop import *
# ReadInput.m

    
    global G_Edof
    command=fscanf(ior,'%s\n',1)
# ReadInput.m:5
    #         disp(command)
    if (strcmp(command,'*ResultDirectory')):
        #             disp('    ResultDirectory')
        command_low=fscanf(ior,'%s\n',1)
# ReadInput.m:9
        Fem.directory = copy(string(command_low))
# ReadInput.m:10
        fprintf(Fem.log,'\t\t*ResultDirectory\n')
        fprintf(Fem.log,'\t\t\t%s\n\n',Fem.directory)
    else:
        if (strcmp(command,'*Nodefile')):
            command_low=fscanf(ior,'%s\n',1)
# ReadInput.m:15
            Node1=Read_Nodefile(command_low)
# ReadInput.m:16
            Node.Coord = copy(Node1(arange(),arange(2,3)))
# ReadInput.m:17
            Node.NNode = copy(size(Node.Coord,1))
# ReadInput.m:18
            Node.Id = copy(int32(Node1(arange(),1)))
# ReadInput.m:19
            Node.F_int = copy(zeros(dot(Node.NNode,2),1))
# ReadInput.m:20
            Node.F_out = copy(zeros(dot(Node.NNode,2),1))
# ReadInput.m:21
            Node.d = copy(zeros(dot(Node.NNode,2),1))
# ReadInput.m:22
            BC_N=zeros(dot(size(Node.Coord,1),2),1)
# ReadInput.m:25
            BC_E=nan(dot(size(Node.Coord,1),2),1)
# ReadInput.m:26
            Node.BC_N = copy(BC_N)
# ReadInput.m:27
            Node.BC_E = copy(BC_E)
# ReadInput.m:28
            fprintf(Fem.log,'\t\t*Nodefile\n')
            fprintf(Fem.log,'\t\t\t%s\n\n',command_low)
        else:
            if (strcmp(command,'*Elementfile')):
                command_low=fscanf(ior,'%s\n',1)
# ReadInput.m:34
                ior1=fopen(command_low,'r')
# ReadInput.m:35
                Element1=Read_Elementfile(command_low)
# ReadInput.m:36
                Element.Connectivity = copy(Element1(arange(),arange(2,end())))
# ReadInput.m:37
                Element.NElem = copy(size(Element.Connectivity,1))
# ReadInput.m:38
                Element.Id = copy(int32(Element1(arange(),1)))
# ReadInput.m:39
                fprintf(Fem.log,'\t\t*Elementfile\n')
                fprintf(Fem.log,'\t\t\t%s\n\n',command_low)
            else:
                if (strcmp(command,'*Elastic')):
                    command_low=fscanf(ior,'%e\t%e\t%e\n',3)
# ReadInput.m:44
                    Fem.E = copy(command_low(1))
# ReadInput.m:45
                    Fem.poisson = copy(command_low(2))
# ReadInput.m:46
                    Fem.thick = copy(command_low(3))
# ReadInput.m:47
                    fprintf(Fem.log,'\t\t*Elastic\n')
                    fprintf(Fem.log,'\t\t\tElastic modulus: %e\n',Fem.E)
                    fprintf(Fem.log,'\t\t\tPoisson ratio: %e\n',Fem.poisson)
                    fprintf(Fem.log,'\t\t\tThickness: %e\n\n',Fem.thick)
                else:
                    if (strcmp(command,'*Density')):
                        command_low=fscanf(ior,'%e\n',1)
# ReadInput.m:54
                        Fem.density = copy(command_low)
# ReadInput.m:55
                        fprintf(Fem.log,'\t\t*Density\n')
                        fprintf(Fem.log,'\t\t\t%e\n\n',command_low)
                    else:
                        if (strcmp(command,'*TimeIntegration')):
                            command=fscanf(ior,'%s\n',1)
# ReadInput.m:60
                            if (strcmp(command,'CDM')):
                                Fem.TimeIntegration = copy('CDM')
# ReadInput.m:62
                            else:
                                if (strcmp(command,'Newmark')):
                                    Fem.TimeIntegration = copy('Newmark')
# ReadInput.m:64
                                    Fem.Newmark[1]=0.25
# ReadInput.m:65
                                    Fem.Newmark[2]=0.5
# ReadInput.m:66
                                    disp('    Newmark coefficient is ')
                                    disp('    beta  = 0.25')
                                    disp('    gamma = 0.5')
                                    disp('    you can change it in ReadInput')
                                else:
                                    if (strcmp(command,'FowardEuler')):
                                        Fem.TimeIntegration = copy('FowardEuler')
# ReadInput.m:72
                                    else:
                                        disp('ForwardEuler, CDM, Newmark are available')
                                        error('Check the Time integration')
                            command_low=fscanf(ior,'%e\t%e\n',2)
# ReadInput.m:78
                            Fem.dt = copy(command_low(1))
# ReadInput.m:79
                            Fem.Time = copy(command_low(2))
# ReadInput.m:80
                            fprintf(Fem.log,'\t\t*TimeIntegration\n')
                            fprintf(Fem.log,'\t\t\t%e\t%e\n\n',command_low(1),command_low(2))
                        else:
                            if (strcmp(command,'*PLE')):
                                command_low=fscanf(ior,'%d\n',1)
# ReadInput.m:85
                                Fem.plane = copy(command_low)
# ReadInput.m:86
                                # Fem.plane == 1 -> plane strain
                                if (Fem.plane == 0):
                                    fprintf(Fem.log,'\t\t*PLE\n')
                                    fprintf(Fem.log,'\t\t\tplane stress condition\n\n')
                                else:
                                    if (Fem.plane == 1):
                                        fprintf(Fem.log,'\t\t*PLE\n')
                                        fprintf(Fem.log,'\t\t\tplane strain condition\n\n')
                                    else:
                                        error('Check *PLE in .inp file')
                                Fem.D = copy(ConstitutiveLaw(Fem))
# ReadInput.m:99
                            else:
                                if (strcmp(command,'*EssentialBD')):
                                    command_low=fscanf(ior,'%s\n',1)
# ReadInput.m:103
                                    ior2=fopen(command_low,'r')
# ReadInput.m:104
                                    if (ior2 == - 1):
                                        error('Could not open the EssentialBD file\n')
                                    reading=textscan(ior2,'%f\t%f\t%f','delimiter',', ')
# ReadInput.m:108
                                    EBDC=concat([reading[1],reading[2],reading[3]])
# ReadInput.m:109
                                    BC_E=nan(dot(size(Node.Coord,1),2),1)
# ReadInput.m:110
                                    BC_E[dot(EBDC(arange(),1),2) + EBDC(arange(),2) - 2]=EBDC(arange(),3)
# ReadInput.m:111
                                    Node.BC_E = copy(BC_E)
# ReadInput.m:112
                                    Node.d = copy(BC_E)
# ReadInput.m:113
                                    fclose(ior2)
                                    fprintf(Fem.log,'\t\t*EssentialBD\n')
                                    fprintf(Fem.log,'\t\t\t%s\n\n',command_low)
                                else:
                                    if (strcmp(command,'*NaturalBD')):
                                        command_low=fscanf(ior,'%s\n',1)
# ReadInput.m:119
                                        ior2=fopen(command_low,'r')
# ReadInput.m:120
                                        if (ior2 == - 1):
                                            error('Could not open the NaturalBD file\n')
                                        reading=textscan(ior2,'%f\t%f\t%f','delimiter',', ')
# ReadInput.m:124
                                        NBDC=concat([reading[1],reading[2],reading[3]])
# ReadInput.m:125
                                        BC_N=zeros(dot(size(Node.Coord,1),2),1)
# ReadInput.m:126
                                        BC_N[dot(NBDC(arange(),1),2) + NBDC(arange(),2) - 2]=NBDC(arange(),3)
# ReadInput.m:127
                                        Node.BC_N = copy(BC_N)
# ReadInput.m:128
                                        fclose(ior2)
                                        fprintf(Fem.log,'\t\t*EssentialBD\n')
                                        fprintf(Fem.log,'\t\t\t%s\n\n',command_low)
                                    else:
                                        if (strcmp(command,'*Step')):
                                            command_low=fscanf(ior,'%s\n',1)
# ReadInput.m:134
                                            fprintf(Fem.log,'\t\t*NaturalBD\n')
                                            fprintf(Fem.log,'\t\t\tto be fixed\n\n')
    