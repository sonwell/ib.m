function [I, Da, Db, DaDa, DaDb, DbDb] = RBFOperators(eps, alpha, beta, gamma, delta)
    psi = @(r, e) 1 ./ sqrt(1 + (e * r).^2);
    dpsi = @(r, e, d) -d .* e^2 ./ (1 + (e * r).^2).^(3/2);
    ddpsi = @(r, e, d1, d2) dpsi(r, e, d1) + 3 * d2 .* e^4 ./ (1 + (e * r).^2).^(5/2);
    phi = @(r) psi(r, eps);
    dphi = @(r, d) dpsi(r, eps, d);
    ddphi = @(r, d1, d2) ddpsi(r, eps, d1, d2);

    ca = DifferenceMatrix(alpha, alpha);
    [ctj, ctk] = ndgrid(beta, beta);
    CM = phi(sqrt(2*(1-sin(ctj).*sin(ctk).*cos(ca)-cos(ctk).*cos(ctj))));

    ea = DifferenceMatrix(gamma, alpha);
    [etj, etk] = ndgrid(delta, beta);
    EDM = sqrt(2*(1-sin(etj).*sin(etk).*cos(ea)-cos(etk).*cos(etj)));

    I = (CM' \ phi(EDM)')';

    Da = (CM' \ dphi(EDM, sin(etj).*sin(etk).*sin(ea))')';
    Db = (CM' \ dphi(EDM, sin(etj).*cos(etk)-cos(etj).*sin(etk).*cos(ea))')';

    DaDa = (CM' \ ddphi(EDM, sin(etj).*sin(etk).*cos(ea), (sin(etj).*sin(etk).*sin(ea)).^2)')';
    DaDb = (CM' \ ddphi(EDM, cos(etj).*sin(etk).*sin(ea), (sin(etj).*sin(etk).*sin(ea)).*(sin(etj).*cos(etk)-cos(etj).*sin(etk).*cos(ea)))')';
    DbDb = (CM' \ ddphi(EDM, cos(etj).*cos(etk)+sin(etj).*sin(etk).*cos(ea), (sin(etj).*cos(etk)-cos(etj).*sin(etk).*cos(ea)).^2)')';
end
