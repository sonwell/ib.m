function [V, l, curr] = Volume(a, b, x, y, z)
    e = 0.892817451684828;
    psi = @(r, e) 1 ./ sqrt(1 + (e * r).^2);
    dpsi = @(r, e, d) -d .* e^2 ./ (1 + (e * r).^2).^(3/2);
    ddpsi = @(r, e, d1, d2) dpsi(r, e, d1) + 3 * d2 .* e^4 ./ (1 + (e * r).^2).^(5/2);
    phi = @(r) psi(r, e);
    dphi = @(r, d) dpsi(r, e, d);
    ddphi = @(r, d1, d2) ddpsi(r, e, d1, d2);

    N = size(a, 1);
    cx = sin(b) .* cos(a);
    cy = sin(b) .* sin(a);
    cz = cos(b);

    ca = DifferenceMatrix(a, a);
    [ctj, ctk] = ndgrid(b, b);
    CDM = sqrt(2*(1-sin(ctj).*sin(ctk).*cos(ca)-cos(ctk).*cos(ctj)));
    CM = phi(CDM);

    EM = CM;
    EDM = CDM;
    etj = ctj;
    etk = ctk;
    ea = ca;

    dEMa = dphi(EDM, sin(etj).*sin(etk).*sin(ea));
    dEMb = dphi(EDM, sin(etj).*cos(etk)-cos(etj).*sin(etk).*cos(ea));

    ddEMaa = ddphi(EDM, sin(etj).*sin(etk).*cos(ea), (sin(etj).*sin(etk).*sin(ea)).^2);
    ddEMab = ddphi(EDM, cos(etj).*sin(etk).*sin(ea), (sin(etj).*sin(etk).*sin(ea)).*(sin(etj).*cos(etk)-cos(etj).*sin(etk).*cos(ea)));
    ddEMbb = ddphi(EDM, cos(etj).*cos(etk)+sin(etj).*sin(etk).*cos(ea), (sin(etj).*cos(etk)-cos(etj).*sin(etk).*cos(ea)).^2);

    pc = property_calculator(CM, dEMa, dEMb, ddEMaa, ddEMab, ddEMbb);
    sphr = pc(cx, cy, cz);
    curr = pc(x, y, z);

    l = sqrt((curr.E .* curr.G - curr.F .^ 2) ./ (sphr.E .* sphr.G - sphr.F .^ 2));
    V = -4 * pi / (3 * N) * l .* dot(curr.x, curr.n, 2);
end

function fn = property_calculator(I, dIa, dIb, dIaa, dIab, dIbb)
    function S = calculate_properties(x, y, z)
        c = I \ [x y z];

        ta = dIa * c;
        tb = dIb * c;

        taa = dIaa * c;
        tab = dIab * c;
        tbb = dIbb * c;

        m = cross(ta, tb, 2);

        E = dot(ta, ta, 2);
        F = dot(ta, tb, 2);
        G = dot(tb, tb, 2);

        J = sqrt(E .* G - F .^ 2);
        n = [1./J 1./J 1./J] .* m;

        e = dot(taa, n, 2);
        f = dot(tab, n, 2);
        g = dot(tbb, n, 2);

        H = (e .* G - 2 * f .* F + g .* E) ./ (2 * (E .* G - F .^ 2));
        K = (e .* g - f .^ 2) ./ (E .* G - F .^ 2);

        S.x = [x y z];
        S.c = c;

        S.ta = ta;
        S.tb = tb;
        S.n = n;

        S.taa = taa;
        S.tab = tab;
        S.tbb = tbb;

        S.E = E; S.e = e;
        S.F = F; S.f = f;
        S.G = G; S.g = g;
        S.H = H;
        S.K = K;

        S.Ea = 2 * dot(taa, ta, 2);
        S.Eb = 2 * dot(tab, ta, 2);
        S.Fa = dot(taa, tb, 2) + dot(ta, tab, 2);
        S.Fb = dot(tab, tb, 2) + dot(ta, tbb, 2);
        S.Ga = 2 * dot(tab, tb, 2);
        S.Gb = 2 * dot(tbb, tb, 2);
    end
    fn = @calculate_properties;
end
