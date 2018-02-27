function [I, Da, Db, Daa, Dab, Dbb, Phi] = RBFOperators(eps, alpha, beta, gamma, delta)
    N = 20;
    meps = 2^-52;  % machine epsilon
    phi = @(r) real(log(r + meps) .* r.^N);
    dphi = @(r, drdx) real((1 + N * log(r + meps)).* r.^(N-2) .* drdx);
    ddphi = @(r, d2rdx2, drdx) real((2 * N - 2 + N * (N-2) * log(r + meps)) .* r.^(N-4) .* drdx + dphi(r, d2rdx2));

    npts = size(alpha, 1);
    ca = DifferenceMatrix(alpha, alpha);
    [ctj, ctk] = ndgrid(beta, beta);
    IDM = sqrt(2*(1-sin(ctj).*sin(ctk).*cos(ca)-cos(ctk).*cos(ctj)));
    ICM = phi(IDM);
    on = ones(npts, 1);
    IM = [ICM on; on' 0];

    mpts = size(gamma, 1);
    ea = DifferenceMatrix(gamma, alpha);
    [etj, etk] = ndgrid(delta, beta);
    EDM = sqrt(2*(1-sin(etj).*sin(etk).*cos(ea)-cos(etk).*cos(etj)));
    om = ones(mpts, 1);
    zm = zeros(mpts, 1);
    EM = [phi(EDM) om];

    IP = (IM' \ EM')';
    I = IP(:, 1:npts);

    DaP = (IM' \ [dphi(EDM, sin(etj).*sin(etk).*sin(ea)) zm]')';
    DbP = (IM' \ [dphi(EDM, sin(etj).*cos(etk)-cos(etj).*sin(etk).*cos(ea)) zm]')';
    Da = DaP(:, 1:npts);
    Db = DbP(:, 1:npts);

    DaaP = (IM' \ [ddphi(EDM, sin(etj).*sin(etk).*cos(ea), (sin(etj).*sin(etk).*sin(ea)).^2) zm]')';
    DabP = (IM' \ [ddphi(EDM, cos(etj).*sin(etk).*sin(ea), (sin(etj).*sin(etk).*sin(ea)).*(sin(etj).*cos(etk)-cos(etj).*sin(etk).*cos(ea))) zm]')';
    DbbP = (IM' \ [ddphi(EDM, cos(etj).*cos(etk)+sin(etj).*sin(etk).*cos(ea), (sin(etj).*cos(etk)-cos(etj).*sin(etk).*cos(ea)).^2) zm]')';
    Daa = DaaP(:, 1:npts);
    Dab = DabP(:, 1:npts);
    Dbb = DbbP(:, 1:npts);

    Phi = IM;
end
