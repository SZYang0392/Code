%Calculate energy density E/V

clear;
clc;
format long;
close all;

%=======================Particle Energy=======================%
%Number of files
fend = 147;
range = 1:fend;

u = 1.660538921e-27;
me = 9.10938291e-31;

% fefi.input
n0 = 1.250228e20;
wtion = [625.1141608304527,        625.1141608304527,     5.420884475596095];
qion = [1, 1, 2];
mi = [2, 3, 4]*u;

wtion = wtion/sum(wtion.*qion, 'all');

%Calculate Energy
Ee = zeros(1, fend);
Ei = zeros(numel(mi), fend);
for k = range
    fpname = ['particle', num2str(k, '%04d'), '.dat'];
    P = loadparticle(fpname);
    Ee(k) = 0.5*me*mean(P.ele.v.^2, 'all');
    for l = 1:numel(mi)
        ind = P.ion.kind == l;
        Ei(l, k) = 0.5*mi(l)*mean(P.ion.v(:, ind).^2, 'all');
    end
end
Ee = Ee*n0;
Ei = Ei*n0.*(wtion.');

%Plot
Epfig = figure;
yyaxis left;
hold on;
plot(range, Ee, 'k.', 'LineWidth', 1.5, 'DisplayName', 'Electron(L)');
Iname = {'Deuterium(L)', 'Tritium(L)', 'Alpha(R)'};
linespec = {'r.', 'm.', 'b.'};
for k = 1:numel(mi)-1
    plot(range, Ei(k, :), linespec{k}, 'LineWidth', 1.5, 'DisplayName', Iname{k});
end
legend('Location', 'NorthWest');
xlabel('FileNum');
ylabel('J/m^3');
k = k + 1;
yyaxis right;
plot(range, Ei(k, :), linespec{k}, 'LineWidth', 1.5, 'DisplayName', Iname{k});
ylabel('J/m^3');
title('Particle Energy');
savefig(Epfig, 'Eparticle.fig');
saveas(Epfig, 'Eparticle.png');


%=======================Field Energy=======================%
epsi0 = 8.854187817e-12;
mu0 = 1.256637061e-6;

%Fefi.input
B0 = [6, 0, 0].'; % in T, [x, y, z];
%Number of files
fend = 147;
range = 1:fend;

Ebw = zeros(1, fend);
Eew = zeros(1, fend);
Efw = zeros(1, fend);

for j = range
    ffname = ['field', num2str(j, '%05d'), '.dat'];
    F = loadfield(ffname);
    Ebw(j) = mean(F.Bx.^2 + F.By.^2 + F.Bz.^2, 'all');
    Eew(j) = mean(F.Ex.^2 + F.Ey.^2 + F.Ez.^2, 'all');
end
% Ebw = (Ebw - sum(B0.^2, 'all'))*0.5/mu0;
Ebw = Ebw*0.5/mu0;
Eew = Eew*0.5*epsi0;

%Plot
Effig = figure;
hold on;
yyaxis left;
plot(range, Eew, 'b.', 'LineWidth', 1.5, 'DisplayName', 'E(L)');
plot(range, Ebw + Eew, 'k.', 'LineWidth', 1.5, 'DisplayName', 'Wave Field(L)');
xlabel('FileNum');
ylabel('J/m^3');
yyaxis right;
plot(range, Ebw, 'r.', 'LineWidth', 1.5, 'DisplayName', 'B(R)');
ylabel('J/m^3');
title('Field Energy');
legend('Location', 'SouthEast');
savefig(Effig, 'Efield.fig');
saveas(Effig, 'Efield.png');