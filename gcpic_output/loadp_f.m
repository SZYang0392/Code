function P = loadp_f(filename)
% 读取文件"particle****.dat"中的粒子数据
% 假设程序的归一化参数采用Gauss制
% 则本程序输出的各个粒子数据都是SI单位制。参数数组除外。

    if isa(filename, 'double')                 %Check data type
        filenamedat= ['particle', num2str(ceil(filename), '%04d'), '.dat'];
        filenamenc = ['particle', num2str(ceil(filename), '%04d'), '.nc'];
    else
        filenamedat = join([filename(1:end-3), 'dat']);
        filenamenc = join([filename(1:end-2), 'nc']);
    end

    if exist(filenamedat, 'file')
        P = loadp_f_dat(filenamedat);
    elseif exist(filenamenc, 'file')
        P = loadp_f_nc(filenamenc);
    else
        fprintf(['File ', filenamedat, ' or ', filenamenc, ' donot exist\n']);
        P.OK = 0;
    end
end

%% Load .dat file
function P = loadp_f_dat(filename)
    % ===================== 打开文件 ===================== %
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
    P.param.stime = pdata(ind);                       %stime，实际上是wci*t
    ind = ind + 1; %3
    P.param.nion = pdata(ind);                        %一共输出了这么多离子，可能不准，后面会修正
    ind = ind + 1; %4
    P.param.nele = pdata(ind);                        %一共输出了这么多电子，可能不准，后面会修正
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
    IndCrit = find(pdata == 125125125);
    if numel(IndCrit) ~= 1
        fprintf(['------------------------ Output Error in ', filename, '------------------------\n']);
    end

    % =====================读取离子数据===================== %
    idata = pdata(ind:IndCrit - 1);
    idata = reshape(idata, P.param.psize, []);
    Ionkind = idata(1,:);                             % 离子种类序号
    P.param.nion = numel(Ionkind);
    % 对离子进行分类处理
    Ionspec0 = tabulate(Ionkind);
    Ionspec = Ionspec0(:, 1);
    if numel(Ionspec) <= P.param.maxkinds
        for k = 1:numel(Ionspec)
            P.ion(k).kind = Ionspec(k);
            P.ion(k).frac = P.param.fraci(k);
            IndI = Ionkind == Ionspec(k);
            P.ion(k).num = sum(IndI);
            P.ion(k).r = idata([2,3,4], IndI).*(P.param.rstep.');
            P.ion(k).v = idata([5,6,7], IndI).*(P.param.vnorm.');
            P.ion(k).mark = idata(8, IndI);
            P.ion(k).cpu = idata(9, IndI);
			P.ion(k).w = idata(end, IndI);
        end
        clear('IndI');
    else
        P.ion.kind = idata(1,:);
        P.ion.frac = idata(11,:);
        P.param.num = size(idata, 2);
        P.ion.r = idata([2,3,4],:).*(P.param.rstep.');
        P.ion.v = idata([5,6,7],:).*(P.param.vnorm.');
        P.ion.mark = idata(8, :);
        P.ion.cpu = idata(9, :);
		P.ion.w = idata(end, :);
    end
    clear('idata', 'Ionkind', 'Ionspec0', 'Ionspec');

    % =====================读取电子数据===================== %
    ind = IndCrit + 1;
    edata = pdata(ind:end);
    edata = reshape(edata, P.param.psize, []);
    Elekind = edata(1, :);
    P.param.nele = numel(Elekind);
    % 分类处理电子
    Elespec0 = tabulate(Elekind);
    Elespec = Elespec0(:, 1);
    if numel(Elespec) <= P.param.maxkinds
        for k = 1:numel(Elespec)
            P.ele(k).kind = Elespec(k);
            P.ele(k).frac = P.param.frace(k);
            IndE = Elekind == Elespec(k);
            P.ele(k).num = sum(IndE);
            P.ele(k).r = edata([2,3,4], IndE).*(P.param.rstep.');
            P.ele(k).v = edata([5,6,7], IndE).*(P.param.vnorm.');
            P.ele(k).mark = edata(8, IndE);
            P.ele(k).cpu = edata(9, IndE);
			P.ele(k).w = edata(end, IndE);
        end
        clear('IndE');
    else
        P.ele.kind = edata(1,:);
        P.ele.frac = edata(11,:);
        P.ele.num = size(edata, 2);
        P.ele.r = edata([2,3,4], :).*(P.param.rstep.');
        P.ele.v = edata([5,6,7],:).*(P.param.vnorm.');
        P.ele.mark = edata(8, :);
        P.ele.cpu = edata(9, :);
		P.ele.w = edata(end, :);
    end
    clear('edata', 'Elekind', 'Elespec0', 'Elespec');
end

%% Load .nc file
function P = loadp_f_nc(filename)
    % ===================== 打开文件 ===================== %
    ncid = netcdf.open(filename,'NC_NOWRITE');
    [~, numvars, ~, ~] = netcdf.inq(ncid);
    data=[];
    for i=1:numvars
        temp = netcdf.getVar(ncid,i-1);
        [name,~,~,~]=netcdf.inqVar(ncid,i-1);
        data.(name) = temp;
    end
    netcdf.close(ncid);
    data.param = data.param.';
    P.OK = 1;
    P.filename = filename;

    % ===================== 整理.nc文件中的数据 ===================== %
    P.nparam =numel(data.param);
    P.paramarray = data.param;
    ind = 1;
    P.param.stime = P.paramarray(ind);                       %stime，实际上是wci*t
    ind = ind + 1; %3
    P.param.nion = P.paramarray(ind);                        %一共输出了这么多离子，可能不准，后面会修正
    ind = ind + 1; %4
    P.param.nele = P.paramarray(ind);                        %一共输出了这么多电子，可能不准，后面会修正
    ind = ind + 1; %5
    P.param.psize = P.paramarray(ind);                       %每个粒子用这么多double存储
    ind = ind + 1; %6
    P.param.maxkinds = P.paramarray(ind);                    %最多多少种电子/离子
    ind = ind + 1; %7
    P.param.betae = P.paramarray(ind);                       %betae
    ind = ind + 1; %8
    P.param.nx = P.paramarray(ind);                          %nx
    ind = ind + 1; %9
    P.param.ny = P.paramarray(ind);                          %ny
    ind = ind + 1; %10
    P.param.nz = P.paramarray(ind);                          %nz
    ind = ind + 1; %11
    P.param.fraci = P.paramarray(ind:ind+P.param.maxkinds-1);%fraci
    ind = ind + 10;%21
    P.param.frace = P.paramarray(ind:ind+P.param.maxkinds-1);%frace
    ind = ind + 10;%31
    P.param.nregion = P.paramarray(ind);                     %nregion
    ind = ind + 1; %32
    xmin = P.paramarray(ind:ind+P.param.nregion-1);
    ind = ind + 10;%42
    xmax = P.paramarray(ind:ind+P.param.nregion-1);
    ind = ind + 10;%52
    ymin = P.paramarray(ind:ind+P.param.nregion-1);
    ind = ind + 10;%62
    ymax = P.paramarray(ind:ind+P.param.nregion-1);
    ind = ind + 10;%72
    zmin = P.paramarray(ind:ind+P.param.nregion-1);
    ind = ind + 10;%82
    zmax = P.paramarray(ind:ind+P.param.nregion-1);
    P.param.pdomain = [xmin; xmax; ymin; ymax; zmin; zmax].';
    % 区间数组，nregion*6，每一行[xmin, xmax, ymin, ymax, zmin, zmax]代表一个区间
    ind = ind + 10;%92
    bx = P.paramarray(ind:ind+P.param.nregion-1);
    ind = ind + 10;%102
    by = P.paramarray(ind:ind+P.param.nregion-1);
    ind = ind + 10;%112
    bz = P.paramarray(ind:ind+P.param.nregion-1);
    P.param.b = [bx, by, bz].';                       % 磁场方向。大小不确定
    % 每一行表示一个区间内的平均磁场方向，列数=观测区间数
    ind = ind + 10;%122
    P.param.rstep = P.paramarray(ind:ind+2)/100.0;
    % in (m)，x、y、z方向一个格点代表多少长度
    ind = ind + 3; %125
    P.param.vnorm = P.paramarray(ind)/100.0;                 % 速度的归一化参数，in (m/s)
    ind = ind + 1; %126
    P.param.va = P.paramarray(ind)/100.0;                    % 归一化Alfven速度？？？？？？？
    
    % =====================读取离子数据===================== %
    if isempty(data.ion)
        P.ion.kind = [];
        P.ion.frac = [];
        P.ion.num = 0;
        P.ion.r = zeros(3, 0);
        P.ion.v = zeros(3, 0);
        P.ion.vb = zeros(3, 0);
        P.ion.w = zeros(1, 0);
        P.ion.ft0 = zeros(1, 0);
        P.ion.mark = zeros(1, 0);
        P.ion.cpu = zeros(1, 0);
        P.ion.marked = zeros(1, 0);
    else
        Ionkind = data.ion(1, :);
        P.param.nion = numel(Ionkind);
        % 对离子进行分类处理
        Ionspec0 = tabulate(Ionkind);
        Ionspec = Ionspec0(:, 1);
        if numel(Ionspec) <= P.param.maxkinds
            for k = 1:numel(Ionspec)
                P.ion(k).kind = Ionspec(k);
                P.ion(k).frac = P.param.fraci(k);
                IndI = Ionkind == Ionspec(k);
                P.ion(k).num = sum(IndI);
                P.ion(k).r = data.ion([2,3,4], IndI).*(P.param.rstep.');
                P.ion(k).v = data.ion([5,6,7], IndI).*(P.param.vnorm.');
                if size(data.ion, 1) >= 15
                    P.ion(k).vb = data.ion([8,9,10], IndI).*(P.param.vnorm.');
                    P.ion(k).w = data.ion(11, IndI);
                    P.ion(k).ft0 = data.ion(12, IndI);
                    P.ion(k).mark = data.ion(13, IndI);
                    P.ion(k).cpu = data.ion(14, IndI);
                    P.ion(k).marked = data.ion(15, IndI);
                else
                    P.ion(k).mark = data.ion(8, IndI);
                    P.ion(k).cpu = data.ion(9, IndI);
                    P.ion(k).w = data.ion(end, IndI);
                end
            end
            clear('IndI');
        else
            P.ion.kind = data.ion(1,:);
            P.ion.frac = data.ion(11,:);
            P.param.num = size(data.ion, 2);
            P.ion.r = data.ion([2,3,4],:).*(P.param.rstep.');
            P.ion.v = data.ion([5,6,7],:).*(P.param.vnorm.');
            if size(data.ion, 1) >= 15
                P.ion.vb = data.ion([8,9,10],:).*(P.param.vnorm.');
                P.ion.w = data.ion(11, :);
                P.ion.ft0 = data.ion(12, :);
                P.ion.mark = data.ion(13, :);
                P.ion.cpu = data.ion(14, :);
                P.ion.marked = data.ion(15, :);
            else
                P.ion.mark = data.ion(8, :);
                P.ion.cpu = data.ion(9, :);
                P.ion.w = data.ion(end, :);
            end
        end
        clear('Ionkind', 'Ionspec0', 'Ionspec');
    end

    % =====================读取电子数据===================== %
    if isempty(data.ele)
        P.ele.kind = [];
        P.ele.frac = [];
        P.ele.num = 0;
        P.ele.r = zeros(3, 0);
        P.ele.v = zeros(3, 0);
        P.ele.vb = zeros(3, 0);
        P.ele.w = zeros(1, 0);
        P.ele.ft0 = zeros(1, 0);
        P.ele.mark = zeros(1, 0);
        P.ele.cpu = zeros(1, 0);
        P.ele.marked = zeros(1, 0);
    else
        Elekind = data.ele(1, :);
        P.param.nele = numel(Elekind);
        % 分类处理电子
        Elespec0 = tabulate(Elekind);
        Elespec = Elespec0(:, 1);
        if numel(Elespec) <= P.param.maxkinds
            for k = 1:numel(Elespec)
                P.ele(k).kind = Elespec(k);
                P.ele(k).frac = P.param.frace(k);
                IndE = Elekind == Elespec(k);
                P.ele(k).num = sum(IndE);
                P.ele(k).r = data.ele([2,3,4], IndE).*(P.param.rstep.');
                P.ele(k).v = data.ele([5,6,7], IndE).*(P.param.vnorm.');
                if size(data.ele, 1) >= 15
                    P.ele(k).vb = data.ele([8,9,10], IndE).*(P.param.vnorm.');
                    P.ele(k).w = data.ele(11, IndE);
                    P.ele(k).ft0 = data.ele(12, IndE);
                    P.ele(k).mark = data.ele(13, IndE);
                    P.ele(k).cpu = data.ele(14, IndE);
                    P.ele(k).marked = data.ele(15, IndE);
                else
                    P.ele(k).mark = data.ele(8, IndE);
                    P.ele(k).cpu = data.ele(9, IndE);
                    P.ele(k).w = data.ele(end, IndE);
                end
            end
            clear('IndE');
        else
            P.ele.kind = data.ele(1,:);
            P.ele.frac = data.ele(11,:);
            P.ele.num = size(data.ele, 2);
            P.ele.r = data.ele([2,3,4], :).*(P.param.rstep.');
            P.ele.v = data.ele([5,6,7],:).*(P.param.vnorm.');
            if size(data.ele, 1) >= 15
                P.ele.vb = data.ele([8,9,10],:).*(P.param.vnorm.');
                P.ele.w = data.ele(11, :);
                P.ele.ft0 = data.ele(12, :);
                P.ele.mark = data.ele(13, :);
                P.ele.cpu = data.ele(14, :);
                P.ele.marked = data.ele(15, :);
            else
                P.ele.mark = data.ele(8, :);
                P.ele.cpu = data.ele(9, :);
                P.ele.w = data.ele(end, :);
            end
        end
        clear('Elekind', 'Elespec0', 'Elespec');
    end
end
