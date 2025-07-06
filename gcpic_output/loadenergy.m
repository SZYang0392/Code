function E = loadenergy(filename, sgn)
    if nargin < 2
        sgn = 'simplify';
    end
    if nargin < 1
        filename = 'energy_t.dat';
    end

    fid = fopen(filename, 'rt');
    fscanf(fid, '      t,EB(3),EE(3),Ept(kinds),E(kinds*5,ex,ey,ez,epar,eper)');
    fData = fscanf(fid,'%f',[1,inf]);
    fclose(fid);

    kinds = fData(2);
    energy = reshape(fData(4:end), 9+7*kinds, []);
    clear('fData');

    E.kinds = kinds;
    E.stime = energy(1, :);
    E.B2 = energy([2,3,4], :);                % Bx^2, By^2, Bz^2
    E.E2 = energy([5,6,7], :);                % Ex^2, Ey^2, Ez^2
    E.EB = sum(E.B2, 1);                      % Magnetic field energy
    E.EE = sum(E.E2, 1);                      % Eletric field energy
    E.Eion = energy(8:7+kinds, :);            % Ion kinetic energy
    E.Eele = energy(8+kinds:7+kinds*2, :);    % Electron kinetic energy
    E.Eionf = energy(8+kinds*2, :);           % Ion fluid energy
    E.Eelef = energy(9+kinds*2, :);           % Ion fluid energy
    E.Pix = energy(10+kinds*2:9+kinds*3, :);  % Ion x momentum
    E.Piy = energy(10+kinds*3:9+kinds*4, :);  % Ion y momentum
    E.Piz = energy(10+kinds*4:9+kinds*5, :);  % Ion z momentum
    E.Ei01 = energy(10+kinds*5:9+kinds*6, :); % ni*Ti0(1)*V
    E.Ei02 = energy(10+kinds*6:9+kinds*7, :); % ni*Ti0(2)*V

    if strcmp(sgn, 'simplify')
        f = fieldnames(E);                    % 读取结构体E的所有字段名称
        for k = 2:numel(f)
            field = E.(f{k});        % 读取结构体字段
            if all(field == 0, 'all')         % 如果结构体字段为0
                E.(f{k}) = [];    % 将该字段设置为空数组
            end
        end
        for k = 7:numel(f)                    % 去除部分数组中全0的行
            field = E.(f{k});
            if isempty(field)
                continue;
            end
            for l = kinds:-1:2
                field1 = field(l, :);
                if all(field1 == 0, 'all')
                    field = field(1:l-1, :);
                end
                E.(f{k}) = field;
            end
        end
    end
end