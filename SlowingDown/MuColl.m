function varargout = MuColl(Pf, Ps, v)
    [mus, mupars, mupers, muepsis] = MuColl_s(Pf, Ps, v);
    mu = sum(mus, 3);
    mupar = sum(mupars, 3);
    muper = sum(mupers, 3);
    muepsi = sum(muepsis, 3);

    varargout{1} = mu;
    varargout{2} = mupar;
    varargout{3} = muper;
    varargout{4} = muepsi;
    varargout{5} = mus;
    varargout{6} = mupars;
    varargout{7} = mupers;
    varargout{8} = muepsis;
end