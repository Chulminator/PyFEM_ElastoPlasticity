
    global G_Edof;


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
