function mesh = ReadMesh(title)
    title1 =sprintf('%s.COH',title);
    fid = fopen(title1,'r');
    if (fid ~= -1)
        fgetl(fid)
        reading1 = textscan(fid,'%d %d %d %d','delimiter',',');
        fclose(fid);
        Cohesive = [reading1{1}, reading1{2}, reading1{3}, reading1{4}];
        mesh.Cohesive = Cohesive;
    end
    
    title1 = sprintf('%s.NODE',title);
    fid = fopen(title1,'r');
    fgetl(fid);
    reading1 = textscan(fid,'%f %f %f %f %f %f %f %f','delimiter',',');
    fclose(fid);
    Node = [reading1{1}, reading1{2}, reading1{3}, reading1{4}, reading1{5}, reading1{6}, reading1{7}, reading1{8}];

    str = sprintf('Number of Node\t: %s',num2str(size(Node,1)));
    disp(str);
    
    
    title1 =sprintf('%s.ELEM',title);
    fid = fopen(title1,'r');
    fgetl(fid);
    [reading1] = textscan(fid,'%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d');
    fclose(fid);

    str = sprintf('Number of Element: %s\n',num2str(size(reading1{1},1)));
    disp(str)
    
    tmp1 = reading1{1};
    tmp2 = reading1{2};
    tmp3 = reading1{3};
    tmp4 = reading1{4};
    tmp5 = reading1{5};
    tmp6 = reading1{6};
    tmp7 = reading1{7};
    tmp8 = reading1{8};
    tmp9 = reading1{9};
    tmp10 = reading1{10};
    tmp11 = reading1{11};
    tmp12 = reading1{12};
    tmp13 = reading1{13};
    tmp14 = reading1{14};
    tmp15 = reading1{15};
    tmp16 = reading1{16};
    tmp17 = reading1{17};
    tmp18 = reading1{18};
    tmp19 = reading1{19};
    tmp20 = reading1{20};
    tmp21 = reading1{21};
    tmp22 = reading1{22};
    tmp23 = reading1{23};
    tmp24 = reading1{24};
    tmp25 = reading1{25};
    tmp26 = reading1{26};
    tmp27 = reading1{27};
    tmp28 = reading1{28};
    tmp29 = reading1{29};
    tmp30 = reading1{30};
    tmp31 = reading1{31};
    tmp32 = reading1{32};
    tmp33 = reading1{33};
    
    tmp34 = [tmp1, tmp2, tmp3, tmp4,...
        tmp5, tmp6, tmp7, tmp8, tmp9, tmp10, tmp11,...
        tmp12, tmp13, tmp14, tmp15, tmp16, tmp17,...
        tmp18, tmp19, tmp20, tmp21, tmp22, tmp23,...
        tmp24, tmp25, tmp26, tmp27, tmp28, tmp29,...
        tmp30, tmp31, tmp32, tmp33];
    A= double(tmp34);
    A(A==0.0) = NaN;

    clear tmp*
    A(A==-1) = NaN;
    mesh.Node     = Node;
    mesh.Element  = A;
end

