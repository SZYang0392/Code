function F = loadfield(filename)
    % 读取文件“field****.dat”中的数据
    % 单位统一为SI制
    if isa(filename, 'double')              %Check data type
        filenamedat = ['field', num2str(ceil(filename), '%05d'), '.dat'];
        filenamenc = ['field', num2str(ceil(filename), '%05d'), '.nc'];
        if exist(filenamedat, 'file')
            F = loadfield_dat(filenamedat);
        elseif exist(filenamenc, 'file')
            F = loadfield_nc(filenamenc);
        else
            fprintf(['File ', filenamedat, ' or ', filenamenc, ' donot exist\n']);
            F.OK = 0;
        end
    else
        if exist(filename, 'file') && strcmp(filename(end-2:end), 'dat')
            F = loadfield_dat(filename);
        elseif exist(filename, 'file') && strcmp(filename(end-1:end), 'nc')
            F = loadfield_nc(filename);
        elseif exist(filename, 'file')
            fprintf(join(['File type ', filename, ' is not supported\n']));
        else
            fprintf(join(['File ', filename, ' donot exist\n']));
        end
    end
end

%% Load .dat file
function F = loadfield_dat(filename)
    % 读取文件“field****.dat”中的数据
    F=[];
    fid = fopen(filename, 'rt');
    if(fid==-1)
        F.OK = 0;                           % 0表示读取失败
        return;
    end
    F.OK = 1;
    F.filename = filename;
    
    % 从文件中读取float放入尺寸为[1, Inf]行向量，Inf表示读到文件末尾
    fData=fscanf(fid,'%f',[1,inf]);
    fclose(fid);

    F.nparam = fData(2);                   % 参数数目
    nxd = fData(3); F.nx = nxd;
    nyd = fData(4); F.ny = nyd;
    nzd = fData(5); F.nz = nzd;
    F.parray = fData(6:6+F.nparam-1);       % diagnosis.f90中的param
    F.stime = F.parray(1);                  % wci*t
    index = 6+F.nparam;
    cord = fData(index:index+nxd+nyd+nzd);
    F.x = (cord(1:nxd)./nxd*(F.parray(35) - F.parray(34)) + F.parray(34))*F.parray(40)*1e-2;    % 坐标网格，单位m
    F.x = reshape(F.x, nxd, 1, 1);
    F.y = (cord(nxd+1:nxd+nyd)./nyd*(F.parray(37) - F.parray(36)) + F.parray(36))*F.parray(40)*1e-2;
    F.y = reshape(F.y, 1, nyd, 1);
    F.z = (cord(nxd+nyd+1:nxd+nyd+nzd)./nzd*(F.parray(39) - F.parray(38)) + F.parray(38))*F.parray(40)*1e-2;
    F.z = reshape(F.z, 1, 1, nzd);
    index = index + nxd + nyd + nzd;
    F.nfield = F.parray(80);                % 总共输出了多少种场

    flength = nxd*nyd*nzd;
    field = reshape(fData(index:end),flength,F.nfield).';
    F.unknown = [];
    unknown = 0;
    for j = 1:F.nfield
        fieldnum = F.parray(75-F.nfield+j);
        dataf = field(j,:);
        switch fieldnum
        case 1
            F.Bx = dataf*F.parray(12)*1e-4;
            F.Bx = reshape(F.Bx, nxd, nyd, nzd);
        case 2
            F.By = dataf*F.parray(12)*1e-4;
            F.By = reshape(F.By, nxd, nyd, nzd);
        case 3
            F.Bz = dataf*F.parray(12)*1e-4;
            F.Bz = reshape(F.Bz, nxd, nyd, nzd);
        case 4
            F.ni = dataf*F.parray(14)*1e6; % in m^-3
            F.ni = reshape(F.ni, nxd, nyd, nzd);
        case 5
            F.ne = dataf*F.parray(14)*1e6;
            F.ne = reshape(F.ne, nxd, nyd, nzd);
        case 6
            F.phi = dataf*F.parray(19)*300; % 电势？
            F.phi = reshape(F.phi, nxd, nyd, nzd);
        case 7
            F.a1x = dataf*F.parray(28)*1e-6; % 磁矢势？
            F.a1x = reshape(F.a1x, nxd, nyd, nzd);
        case 8
            F.a1y = dataf*F.parray(28)*1e-6; % 磁矢势？
            F.a1y = reshape(F.a1y, nxd, nyd, nzd);
        case 9
            F.a1z = dataf*F.parray(28)*1e-6; % 磁矢势？
            F.a1z = reshape(F.a1z, nxd, nyd, nzd);
        case 12
            F.Ex = dataf*F.parray(13);
            F.Ex = reshape(F.Ex, nxd, nyd, nzd);
        case 13
            F.Ey = dataf*F.parray(13);
            F.Ey = reshape(F.Ey, nxd, nyd, nzd);
        case 14
            F.Ez = dataf*F.parray(13);
            F.Ez = reshape(F.Ez, nxd, nyd, nzd);
        case 15
            F.uix = dataf*F.parray(15)*1e-2; % 流体速度
            F.uix = reshape(F.uix, nxd, nyd, nzd);
        case 16
            F.uiy = dataf*F.parray(15)*1e-2;%----------------
            F.uiy = reshape(F.uiy, nxd, nyd, nzd);
        case 17
            F.uiz = dataf*F.parray(15)*1e-2;%----------------
            F.uiz = reshape(F.uiz, nxd, nyd, nzd);
        case 18
            F.uex = dataf*F.parray(15)*1e-2;%----------------
            F.uex = reshape(F.uex, nxd, nyd, nzd);
        case 19
            F.uey = dataf*F.parray(15)*1e-2;%----------------
            F.uey = reshape(F.uey, nxd, nyd, nzd);
        case 20
            F.uez = dataf*F.parray(15)*1e-2;%----------------
            F.uez = reshape(F.uez, nxd, nyd, nzd);
        case 21
            F.dJx = dataf*F.parray(18)/(3e5);%----------------
            F.dJx = reshape(F.dJx, nxd, nyd, nzd);
        case 22
            F.dJy = dataf*F.parray(18)/(3e5);%----------------
            F.dJy = reshape(F.dJy, nzd, nyd, nzd);
        case 23
            F.dJz = dataf*F.parray(18)/(3e5);%----------------
            F.dJz = reshape(F.dJz, nzd, nyd, nzd);
        case 24
            F.dBx = dataf*F.parray(12)*1e-4;%----------------
            F.dBx = reshape(F.dBx, nzd, nyd, nzd);
        case 25
            F.dBy = dataf*F.parray(12)*1e-4;%----------------
            F.dBy = reshape(F.dBy, nzd, nyd, nzd);
        case 26
            F.dBz = dataf*F.parray(12)*1e-4;%----------------
            F.dBz = reshape(F.dBz, nzd, nyd, nzd);
        case 27
            F.dEx = dataf*F.parray(13)*3e4;%----------------
            F.dEx = reshape(F.dEx, nzd, nyd, nzd);
        case 28
            F.dEy = dataf*F.parray(13)*3e4;%----------------
            F.dEy = reshape(F.dEy, nzd, nyd, nzd);
        case 29
            F.dEz = dataf*F.parray(13)*3e4;%----------------
            F.dEz = reshape(F.dEz, nzd, nyd, nzd);
        otherwise
            fprintf('Unrecognized field #%d\n', fieldnum);
            unknown = unknown + 1;
            F.unknown(unknown) = fieldnum;
        end
    end
end

%% Load .nc file
function F = loadfield_nc(filename)
    % 读取文件"field****.nc"中的数据
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
    F.OK = 1;
    F.filename = filename;

    % ===================== 关键参数 ===================== %
    F.nparam = numel(data.param);                   % 参数数目
    nxd = numel(data.X); F.nx = nxd;
    nyd = numel(data.Y); F.ny = nyd;
    nzd = numel(data.Z); F.nz = nzd;
    F.parray = data.param;                          % diagnosis_nc.f90中的param
    F.stime = F.parray(1);                          % wci*t
    F.x = (data.X./nxd*(F.parray(35) - F.parray(34)) + F.parray(34))*F.parray(40)*1e-2;
    F.x = reshape(F.x, [], 1, 1);
    F.y = (data.Y./nyd*(F.parray(37) - F.parray(36)) + F.parray(36))*F.parray(40)*1e-2;
    F.y = reshape(F.y, 1, [], 1);
    F.z = (data.Z./nzd*(F.parray(39) - F.parray(38)) + F.parray(38))*F.parray(40)*1e-2;
    F.z = reshape(F.z, 1, 1, []);
    F.nfield = F.parray(80);                      % 总共输出了多少种场
    F.unknown = [];

    % ===================== 单位转换为国际单位制 ===================== %
    F.filename = filename;
    fields = {'param', 'X', 'Y', 'Z'};
    if isfield(data, 'Bx')
        fields = [fields, {'Bx'}];
        F.Bx = data.Bx*F.parray(12)*1e-4;
        F.Bx = reshape(F.Bx, nxd, nyd, nzd);
    end
    if isfield(data, 'By')
        fields = [fields, {'By'}];
        F.By = data.By*F.parray(12)*1e-4;
        F.By = reshape(F.By, nxd, nyd, nzd);
    end
    if isfield(data, 'Bz')
        fields = [fields, {'Bz'}];
        F.Bz = data.Bz*F.parray(12)*1e-4;
        F.Bz = reshape(F.Bz, nxd, nyd, nzd);
    end
    if isfield(data, 'Ni')
        fields = [fields, {'Ni'}];
        F.ni = data.Ni*F.parray(14)*1e6; % in m^-3
        F.ni = reshape(F.ni, nxd, nyd, nzd);
    end
    if isfield(data, 'Ne')
        fields = [fields, {'Ne'}];
        F.ne = data.Ne*F.parray(14)*1e6; % in m^-3
        F.ne = reshape(F.ne, nxd, nyd, nzd);
    end
    if isfield(data, 'phi')
        fields = [fields, {'phi'}];
        F.phi = data.phi*F.parray(19)*300; % 电势？
        F.phi = reshape(F.phi, nxd, nyd, nzd);
    end
    if isfield(data, 'A1x')
        fields = [fields, {'A1x'}];
        F.a1x = data.A1x*F.parray(28)*1e-6; % 磁矢势？
        F.a1x = reshape(F.a1x, nxd, nyd, nzd);
    end
    if isfield(data, 'A1y')
        fields = [fields, {'A1y'}];
        F.a1y = data.A1y*F.parray(28)*1e-6; % 磁矢势？
        F.a1y = reshape(F.a1y, nxd, nyd, nzd);
    end
    if isfield(data, 'A1z')
        fields = [fields, {'A1z'}];
        F.a1z = data.A1z*F.parray(28)*1e-6; % 磁矢势？
        F.a1z = reshape(F.a1z, nxd, nyd, nzd);
    end
    if isfield(data, 'Ex')
        fields = [fields, {'Ex'}];
        F.Ex = data.Ex*F.parray(13);
        F.Ex = reshape(F.Ex, nxd, nyd, nzd);
    end
    if isfield(data, 'Ey')
        fields = [fields, {'Ey'}];
        F.Ey = data.Ey*F.parray(13);
        F.Ey = reshape(F.Ey, nxd, nyd, nzd);
    end
    if isfield(data, 'Ez')
        fields = [fields, {'Ez'}];
        F.Ez = data.Ez*F.parray(13);
        F.Ez = reshape(F.Ez, nxd, nyd, nzd);
    end
    if isfield(data, 'Upx')
        fields = [fields, {'Upx'}];
        F.uix = data.Upx*F.parray(15)*1e-2; % 流体速度
        F.uix = reshape(F.uix, nxd, nyd, nzd);
    end
    if isfield(data, 'Upy')
        fields = [fields, {'Upy'}];
        F.uiy = data.Upy*F.parray(15)*1e-2;%----------------
        F.uiy = reshape(F.uiy, nxd, nyd, nzd);
    end
    if isfield(data, 'Upz')
        fields = [fields, {'Upz'}];
        F.uiz = data.Upz*F.parray(15)*1e-2;%----------------
        F.uiz = reshape(F.uiz, nxd, nyd, nzd);
    end
    if isfield(data, 'Uex')
        fields = [fields, {'Uex'}];
        F.uex = data.Uex*F.parray(15)*1e-2;%----------------
        F.uex = reshape(F.uex, nxd, nyd, nzd);
    end
    if isfield(data, 'Uey')
        fields = [fields, {'Uey'}];
        F.uey = data.Uey*F.parray(15)*1e-2;%----------------
        F.uey = reshape(F.uey, nxd, nyd, nzd);
    end
    if isfield(data, 'Uez')
        fields = [fields, {'Uez'}];
        F.uez = data.Uez*F.parray(15)*1e-2;%----------------
        F.uez = reshape(F.uez, nxd, nyd, nzd);
    end
    if isfield(data, 'Jx')
        fields = [fields, {'Jx'}];
        F.dJx = data.Jx*F.parray(18)/(3e5);%----------------
        F.dJx = reshape(F.dJx, nxd, nyd, nzd);
    end
    if isfield(data, 'Jy')
        fields = [fields, {'Jy'}];
        F.dJy = data.Jy*F.parray(18)/(3e5);%----------------
        F.dJy = reshape(F.dJy, nzd, nyd, nzd);
    end
    if isfield(data, 'Jz')
        fields = [fields, {'Jz'}];
        F.dJz = data.Jz*F.parray(18)/(3e5);%----------------
        F.dJz = reshape(F.dJz, nzd, nyd, nzd);
    end
    if isfield(data, 'B1x')
        fields = [fields, {'B1x'}];
        F.dBx = data.B1x*F.parray(12)*1e-4;%----------------
        F.dBx = reshape(F.dBx, nzd, nyd, nzd);
    end
    if isfield(data, 'B1y')
        fields = [fields, {'B1y'}];
        F.dBy = data.B1y*F.parray(12)*1e-4;%----------------
        F.dBy = reshape(F.dBy, nzd, nyd, nzd);
    end
    if isfield(data, 'B1z')
        fields = [fields, {'B1z'}];
        F.dBz = data.B1z*F.parray(12)*1e-4;%----------------
        F.dBz = reshape(F.dBz, nzd, nyd, nzd);
    end
    if isfield(data, 'E1x')
        fields = [fields, {'E1x'}];
        F.dEx = data.E1x*F.parray(13)*3e4;%----------------
        F.dEx = reshape(F.dEx, nzd, nyd, nzd);
    end
    if isfield(data, 'E1y')
        fields = [fields, {'E1y'}];
        F.dEy = data.E1y*F.parray(13)*3e4;%----------------
        F.dEy = reshape(F.dEy, nzd, nyd, nzd);
    end
    if isfield(data, 'E1z')
        fields = [fields, {'E1z'}];
        F.dEz = data.E1z*F.parray(13)*3e4;%----------------
        F.dEz = reshape(F.dEz, nzd, nyd, nzd);
    end

    % ===================== 读取未知输出 ===================== %
    S = fieldnames(data);
    for j = 1:numel(S)
        IsField = false;
        for k = 1:numel(fields)
            if strcmp(S{j}, fields{k})
                IsField = true;
                break;
            end
        end
        if ~IsField
            F.unknown.(S{j}) = data.(S{j});
        end
    end
end