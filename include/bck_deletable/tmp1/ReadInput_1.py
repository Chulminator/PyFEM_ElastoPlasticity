# Generated with SMOP  0.41
from libsmop import *
# ReadInput.m

    
    global G_Edof
    disp(input_name)
    disp('#ReadInput')
    Node=struct('Coord',1)
# ReadInput.m:6
    Element=struct('Connectivity',1)
# ReadInput.m:7
    ior=fopen(input_name,'r')
# ReadInput.m:8
    if (ior == - 1):
        disp('Could not open the input file\n')
        error('FE analysis code is terminated, readInput')
    
    command=fscanf(ior,'%s\n',1)
# ReadInput.m:14
    if (strcmp(command,'*Title')):
        command_low=fscanf(ior,'%s\n',1)
# ReadInput.m:16
        Fem=struct('title',command_low)
# ReadInput.m:17
    else:
        error('The first information in .inp should be *Title')
    
    
    tmp_str=sprintf('./log/log_%s.txt',Fem.title)
# ReadInput.m:22
    Fem.log = copy(fopen(tmp_str,'w'))
# ReadInput.m:23
    fprintf(Fem.log,input_name)
    fprintf(Fem.log,'\n\t#ReadInput\n')
    fprintf(Fem.log,'\t\t*Title\n')
    fprintf(Fem.log,'\t\t\t%s\n\n',Fem.title)