function [S, D, P] = SDP(w, B, Ps)
    Sa = zeros(numel(Ps), numel(Ps(1).n0));
    Da = zeros(numel(Ps), numel(Ps(1).n0));
    Pa = zeros(numel(Ps), numel(Ps(1).n0));
    for k = 1:numel(Ps)
        [Sa(k, :), Da(k, :), Pa(k, :)] = KhiEMCold(Ps(k), B, w);
    end
    S = sum(Sa, 1);
    D = sum(Da, 1);
    P = sum(Pa, 1);
end