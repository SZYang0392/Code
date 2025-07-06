function P = loadp(filename)
% 读取文件"particle****.dat"中的粒子数据
% 假设程序的归一化参数采用Gauss制
% 则本程序输出的各个粒子数据都是SI单位制。参数数组除外。

% ！！！！diagnosis_yang.f90里面按照fraci来区分粒子种类，换言之fraci相同的粒子种类无法区分

    if isa(filename, 'double')                 %Check data type
        filename = ['particle', num2str(ceil(filename), '%04d'), '.dat'];
    end

    fid = fopen(filename, 'rt');
    if (fid == -1)
        P.OK = 0;                              % 0表示读取失败
        fprintf(['cannot open ', filename, '\n']);
        return;
    end
    P.OK = 1;                                  % 成功读取文件
    P.filename = filename;

    % 从文件中读取float放入尺寸为[1, Inf]行向量，Inf表示读到文件末尾
    pdata = fscanf(fid, '%f', [1, Inf]);
    fclose(fid);

    ind = 1;                                   %目前该读pdata的第几号数据
    % =====================读取程序中的param数组===================== %
    P.nparam = pdata(ind);                            %参数数组param一共这么大
    P.paramarray = pdata(ind+1 : ind+P.nparam);       %数组param
    ind = ind + 1; %2
    P.param.stime = pdata(ind);                       %stime，大师姐的diagnosis_nc.f90中说归一化到omega_i
    ind = ind + 1; %3
    P.param.nion = pdata(ind);                        %一共输出了这么多离子
    ind = ind + 1; %4
    P.param.nele = pdata(ind);                        %一共输出了这么多电子
    ind = ind + 1; %5
    P.param.psize = pdata(ind);                       %每个粒子用这么多double存储
    ind = ind + 1; %6
    P.param.maxkinds = pdata(ind);                    %最多多少种电子/离子
    ind = ind + 1; %7
    P.param.betae = pdata(ind);                       %betae
    ind = ind + 1; %8
    P.param.nx = pdata(ind);                          %nx
    ind = ind + 1; %9
    P.param.ny = pdata(ind);                          %ny
    ind = ind + 1; %10
    P.param.nz = pdata(ind);                          %nz
    ind = ind + 1; %11
    P.param.fraci = pdata(ind:ind+P.param.maxkinds-1);%fraci
    ind = ind + 10;%21
    P.param.frace = pdata(ind:ind+P.param.maxkinds-1);%frace
    ind = ind + 10;%31
    P.param.nregion = pdata(ind);                     %nregion
    ind = ind + 1; %32
    xmin = pdata(ind:ind+P.param.nregion-1);
    ind = ind + 10;%42
    xmax = pdata(ind:ind+P.param.nregion-1);
    ind = ind + 10;%52
    ymin = pdata(ind:ind+P.param.nregion-1);
    ind = ind + 10;%62
    ymax = pdata(ind:ind+P.param.nregion-1);
    ind = ind + 10;%72
    zmin = pdata(ind:ind+P.param.nregion-1);
    ind = ind + 10;%82
    zmax = pdata(ind:ind+P.param.nregion-1);
    P.param.pdomain = [xmin; xmax; ymin; ymax; zmin; zmax].';
    % 区间数组，nregion*6，每一行[xmin, xmax, ymin, ymax, zmin, zmax]代表一个区间
    ind = ind + 10;%92
    bx = pdata(ind:ind+P.param.nregion-1);
    ind = ind + 10;%102
    by = pdata(ind:ind+P.param.nregion-1);
    ind = ind + 10;%112
    bz = pdata(ind:ind+P.param.nregion-1);
    P.param.b = [bx, by, bz].';                       % 磁场方向。大小不确定
    % 每一行表示一个区间内的平均磁场方向，列数=观测区间数
    ind = ind + 10;%122
    P.param.rstep = pdata(ind:ind+2)/100.0;
    % in (m)，x、y、z方向一个格点代表多少长度
    ind = ind + 3; %125
    P.param.vnorm = pdata(ind)/100.0;                 % 速度的归一化参数，in (m/s)
    ind = ind + 1; %126
    P.param.va = pdata(ind)/100.0;                    % 归一化Alfven速度？？？？？？？
    ind = ind + 1;

    % =====================检查输出是否完整===================== %
    Nok = 2 + P.nparam + (P.param.nele + P.param.nion)*P.param.psize;
    npdata = numel(pdata);
    if npdata < Nok
        fprintf(['---------------------Insufficient output in', filename, ' ---------------------\n']);
    elseif npdata > Nok
        fprintf(['---------------------Extra output in', filename, ' ---------------------\n']);
    end

    % =====================读取离子数据===================== %
    idata = pdata(ind:ind+P.param.nion*P.param.psize-1);
    idata = reshape(idata, P.param.psize, P.param.nion);
    P.ion.kind = idata(1,:);                          % 离子种类序号
    % 离子位置数组，单位是m
    P.ion.r = idata([2,3,4],:);
    for k = 1:3
        P.ion.r(k,:) = P.ion.r(k,:)*P.param.rstep(k);
    end
    P.ion.v = idata([5,6,7],:).*(P.param.vnorm.');
    % 离子速度数组，单位m/s
    P.ion.vb = idata([8,9,10],:).*(P.param.vnorm.');
    % 离子速度数组，单位m/s，但是按照磁场取直角坐标系。第一列是平行磁场方向
    P.ion.frac = idata(11,:);
    % 离子权重数组
    clear('idata');

    % =====================读取电子数据===================== %
    ind = ind + P.param.nion*P.param.psize;
    
    if pdata(ind) == 125125125
        fprintf('Ion OK\n');
    end
    ind = ind + 1

    edata = pdata(ind:ind+P.param.nele*P.param.psize-1);
    edata = reshape(edata, P.param.psize, P.param.nele);
    P.ele.kind = edata(1,:);
    P.ele.r = edata([2,3,4],:);
    for k = 1:3
        P.ele.r(k,:) = P.ele.r(k,:)*P.param.rstep(k);
    end
    % 电子位置数组，单位是m
    P.ele.v = edata([5,6,7],:).*(P.param.vnorm.');
    % 电子速度数组，单位m/s
    P.ele.vb = edata([8,9,10],:).*(P.param.vnorm.');
    % 电子速度数组，单位m/s，但是按照磁场取直角坐标系。第一列是平行磁场方向
    P.ele.frac = edata(11,:);
    % 电子权重数组
    clear('edata');
end