function [Node, Element, Fem] =ReadInput(input_name)
    global G_Edof;

    disp(input_name)
    disp('#ReadInput')
    Node = struct('Coord',1);
    Element= struct('Connectivity',1);
    ior = fopen(input_name,'r');
    if(ior == -1)
        disp('Could not open the input file\n');
        error('FE analysis code is terminated, readInput');
    end

    command = fscanf(ior, '%s\n',1);
    if (strcmp(command,"*Title"))
        command_low = fscanf(ior, '%s\n',1);
        Fem = struct('title',command_low);
    else
        error('The first information in .inp should be *Title');
    end
    
    tmp_str = sprintf('./log/log_%s.txt',Fem.title);
    Fem.log = fopen(tmp_str,'w');
    fprintf(Fem.log,input_name);
    fprintf(Fem.log,'\n\t#ReadInput\n');

    fprintf(Fem.log,'\t\t*Title\n');
    fprintf(Fem.log,'\t\t\t%s\n\n', Fem.title);

    while ~feof(ior)
        command = fscanf(ior,'%s\n',1);
%         disp(command)
        if (strcmp(command,"*ResultDirectory"))
%             disp('    ResultDirectory')
            command_low = fscanf(ior, '%s\n',1);
            Fem.directory=string(command_low);
            fprintf(Fem.log,'\t\t*ResultDirectory\n');
            fprintf(Fem.log,'\t\t\t%s\n\n', Fem.directory);
    
        elseif (strcmp(command,"*Nodefile"))
            command_low = fscanf(ior, '%s\n',1);
            [Node1] =Read_Nodefile(command_low);
            Node.Coord = Node1(:,2:3);
            Node.NNode = size(Node.Coord,1);
            Node.Id = int32(Node1(:,1));
            Node.F_int= zeros(Node.NNode*2,1);
            Node.F_out= zeros(Node.NNode*2,1);
            Node.d= zeros(Node.NNode*2,1);     %% displacement


            BC_N = zeros(size(Node.Coord,1)*2,1);
            BC_E = nan(size(Node.Coord,1)*2,1);
            Node.BC_N = BC_N;
            Node.BC_E = BC_E;

            fprintf(Fem.log,'\t\t*Nodefile\n');
            fprintf(Fem.log,'\t\t\t%s\n\n', command_low);
    
        elseif (strcmp(command,"*Elementfile"))
            command_low = fscanf(ior, '%s\n',1);
            ior1 = fopen(command_low,'r');
            [Element1] =Read_Elementfile(command_low);
            Element.Connectivity = Element1(:,2:end);
            Element.NElem = size(Element.Connectivity, 1);
            Element.Id = int32(Element1(:,1));
            fprintf(Fem.log,'\t\t*Elementfile\n');
            fprintf(Fem.log,'\t\t\t%s\n\n', command_low);
    
        elseif (strcmp(command,"*Elastic"))
            command_low = fscanf(ior,'%e\t%e\t%e\n',3);
            Fem.E=command_low(1);
            Fem.poisson=command_low(2);
            Fem.thick=command_low(3);
            fprintf(Fem.log,'\t\t*Elastic\n');
            fprintf(Fem.log,'\t\t\tElastic modulus: %e\n', Fem.E);
            fprintf(Fem.log,'\t\t\tPoisson ratio: %e\n', Fem.poisson);
            fprintf(Fem.log,'\t\t\tThickness: %e\n\n', Fem.thick);
    
        elseif (strcmp(command,"*Density"))
            command_low = fscanf(ior,'%e\n',1);
            Fem.density = command_low;
            fprintf(Fem.log,'\t\t*Density\n');
            fprintf(Fem.log,'\t\t\t%e\n\n', command_low);
    
        elseif (strcmp(command,"*TimeIntegration"))
            command = fscanf(ior,'%s\n',1);
            if(strcmp(command, "CDM"))
                Fem.TimeIntegration = "CDM";
            elseif(strcmp(command, "Newmark"))
                Fem.TimeIntegration = "Newmark";
                Fem.Newmark(1) = 0.25;
                Fem.Newmark(2) = 0.5;
                disp("    Newmark coefficient is ")
                disp("    beta  = 0.25")
                disp("    gamma = 0.5")
                disp("    you can change it in ReadInput")
            elseif(strcmp(command, "FowardEuler"))
                Fem.TimeIntegration = "FowardEuler";
            else
                disp("ForwardEuler, CDM, Newmark are available")
                error("Check the Time integration")
            end
            
            command_low = fscanf(ior,'%e\t%e\n',2);
            Fem.dt  = command_low(1);
            Fem.Time= command_low(2);
            fprintf(Fem.log,'\t\t*TimeIntegration\n');
            fprintf(Fem.log,'\t\t\t%e\t%e\n\n', command_low(1), command_low(2));
    
        elseif (strcmp(command,"*PLE"))
            command_low = fscanf(ior,'%d\n',1);
            Fem.plane = command_low;
            % Fem.plane == 0 -> plane stress
            % Fem.plane == 1 -> plane strain
            if(Fem.plane == 0)
                fprintf(Fem.log,'\t\t*PLE\n');
                fprintf(Fem.log,'\t\t\tplane stress condition\n\n');
            elseif(Fem.plane == 1)
                fprintf(Fem.log,'\t\t*PLE\n');
                fprintf(Fem.log,'\t\t\tplane strain condition\n\n');
            else
                error('Check *PLE in .inp file');
            end
            
            Fem.D = ConstitutiveLaw(Fem);


        elseif (strcmp(command,"*EssentialBD"))
            command_low = fscanf(ior,'%s\n',1);
            ior2 = fopen(command_low,'r');
            if(ior2 == -1)
                error('Could not open the EssentialBD file\n');
            end
            reading = textscan(ior2,'%f\t%f\t%f','delimiter',', ');
            EBDC = [reading{1}, reading{2}, reading{3}];
            BC_E = nan(size(Node.Coord,1)*2,1);
            BC_E(EBDC(:,1)*2+EBDC(:,2)-2)=EBDC(:,3);
            Node.BC_E = BC_E;
            Node.d = BC_E;
            fclose(ior2);

            fprintf(Fem.log,'\t\t*EssentialBD\n');
            fprintf(Fem.log,'\t\t\t%s\n\n',command_low);
        elseif (strcmp(command,"*NaturalBD"))
            command_low = fscanf(ior,'%s\n',1);
            ior2 = fopen(command_low,'r');
            if(ior2 == -1)
                error('Could not open the NaturalBD file\n');
            end
            reading = textscan(ior2,'%f\t%f\t%f','delimiter',', ');
            NBDC = [reading{1}, reading{2}, reading{3}];
            BC_N = zeros(size(Node.Coord,1)*2,1);
            BC_N(NBDC(:,1)*2+NBDC(:,2)-2)=NBDC(:,3);
            Node.BC_N = BC_N;
            fclose(ior2);
            fprintf(Fem.log,'\t\t*EssentialBD\n');
            fprintf(Fem.log,'\t\t\t%s\n\n',command_low);
    
        elseif (strcmp(command,"*Step"))  %% to be fixed
            command_low = fscanf(ior,'%s\n',1);

            fprintf(Fem.log,'\t\t*NaturalBD\n');
            fprintf(Fem.log,'\t\t\tto be fixed\n\n');
        end
    end

    Fem = DefineGP(Fem);

%     Analysis mode
    if (1)
        addpath("include\Dynamics")
        % dynamic property
        Node.v = zeros(Node.NNode*2,1);     %% Acceleration
        Node.a = zeros(Node.NNode*2,1);     %% displacement
        Node.M = zeros(Node.NNode*2,1);       %% Mass
%         Node.M = zeros(Node.NNode,1);       %% Mass
        Node.Fint= zeros(Node.NNode*2,1);
        Node.Fext= zeros(Node.NNode*2,1);
        Node.Ftotal= zeros(Node.NNode*2,1);
    end
    if (0)
        % UP formulation
        Node.BC_Ep = zeros(Node.NNode,1);  
        Node.BC_Np = zeros(Node.NNode,1);  
        Node.p = zeros(Node.NNode,1);  
    end
    if (1)
        addpath("include\PhaseField")
        Fem.PhaseFieldFlag = 1;
        % Phase field
        Element.H = zeros(Node.NNode,Fem.Ngp);   %% history variable
        Element.Strain = zeros(Node.NNode,Fem.Ngp);   %% Strain at Gauss Point
%         Element.Kelephase = zeros(Node.NNode,1);    %% Stiffness matrix for phase field
        Node.Rextphase = zeros(Node.NNode,1);    %% Rext phase
        Node.Phi = zeros(Node.NNode,1);    %% Damage variabl
        % e
        [V, C] = eig(Fem.D);

        LambdaPlus = C*(C>0); 
        LambdaMinus = C - LambdaPlus;    
    
        DPlus = V * LambdaPlus * inv(V); % strain 
        DMinus = V * LambdaMinus * inv(V);    

        Fem.DPlus = DPlus;
        Fem.DMinus = DMinus;
    end


    fclose(ior);
end

function [Node] =Read_Nodefile(NodefileStr)
    global G_Dim;

    ior = fopen(NodefileStr,'r');
    if(ior==-1)
        error('Check the node file - name or direction');
    end
    NNode = str2num(fgetl(ior));
    if(G_Dim == 2)
        reading = textscan(ior,'%f\t%f\t%f','delimiter',', ');
        Node=[reading{1},reading{2},reading{3}];
        if(NNode ~= size(Node,1))
            disp(NodefileStr)
            error("Check here!")
        end
        fclose(ior);
    elseif(G_Dim ==3)
        error("3D is not implemented YET")
    end
end

function [Element] =Read_Elementfile(ElemfileStr)
    global G_Edof;

    ior = fopen(ElemfileStr,'r');
    if(ior==-1)
        error('Check the Element file - name or direction');
    end
    NElement = str2num(fgetl(ior));
    if(G_Edof == 4)
        reading = textscan(ior,'%d\t%d\t%d\t%d\t%d');
        Element=[reading{1},reading{2},reading{3},reading{4},reading{5}];
        if(NElement ~= size(Element,1))
            disp(ElemfileStr)
            error("Check here!")
        end
        fclose(ior);

    elseif(G_Edof == 3)
        reading = textscan(ior,'%d\t%d\t%d\t%d');
        Element=[reading{1},reading{2},reading{3},reading{4}];
        if(NElement ~= size(Element,1))
            disp(ElemfileStr)
            error("Check here!")
        end
        fclose(ior);

    elseif(G_Edof == 6)
        reading = textscan(ior,'%d\t%d\t%d\t%d\t%d\t%d\t%d');
        Element=[reading{1},reading{2},reading{3},reading{4},reading{5},reading{6},reading{7}];
        if(NElement ~= size(Element,1))
            disp(ElemfileStr)
            error("Check here!")
        end
        fclose(ior);
    else
        error("T3, T6, Q4 eleement is available for now")
    end    
end

function Fem = DefineGP(Fem)
global G_Edof;
    if(G_Edof == 4)
        Fem.Ngp = 4;
        Fem.GP = [-sqrt(1.0/3.0) -sqrt(1.0/3.0) 1.0;
                   sqrt(1.0/3.0) -sqrt(1.0/3.0) 1.0;
                   sqrt(1.0/3.0)  sqrt(1.0/3.0) 1.0;
                  -sqrt(1.0/3.0)  sqrt(1.0/3.0) 1.0];
%         Fem.Weight = [1.0 1.0 1.0 1.0];

    elseif(G_Edof == 3)
        Fem.Ngp = 1;
        Fem.GP = [1.0/3.0 1.0/3.0 0.5];

    elseif(G_Edof == 6)
        Fem.Ngp = 3;
        Fem.GP = [1.0/6.0 1.0/6.0 1.0/6.0
                  1.0/6.0 2.0/3.0 1.0/6.0
                  2.0/3.0 1.0/6.0 1.0/6.0];
    else
        error("T3, T6, Q4 element is available")
    end
end
