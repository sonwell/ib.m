function f = ForceCalculator(alpha, beta, x0, y0, z0, youngs, bulk, bending)
    f = force_calculator(1.892817451684828, alpha, beta, x0, y0, z0, youngs, bulk, bending);
end

function calculator = force_calculator(e, a, b, x0, y0, z0, youngs, bulk, bending)
    psi = @(r, e) 1 ./ sqrt(1 + (e * r).^2);
    dpsi = @(r, e, d) -d .* e^2 ./ (1 + (e * r).^2).^(3/2);
    ddpsi = @(r, e, d1, d2) dpsi(r, e, d1) + 3 * d2 .* e^4 ./ (1 + (e * r).^2).^(5/2);
    phi = @(r) psi(r, e);
    dphi = @(r, d) dpsi(r, e, d);
    ddphi = @(r, d1, d2) ddpsi(r, e, d1, d2);


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
    orig = pc(x0, y0, z0);
    
    function [fsk, fb, sdata, bdata]= f_calculator(x, y, z)
        curr = pc(x, y, z);
        [fsk, sdata] = sk_force(orig, curr, youngs, bulk);
        [fb, bdata] = bending_force(orig, curr, bending);
    end
    calculator = @f_calculator;
end

function [f, data] = sk_force(orig, curr, E, b)
    detG = curr.E .* curr.G - curr.F .^ 2;
    detGa = curr.Ea .* curr.G + curr.E .* curr.Ga - 2 * curr.F .* curr.Fa;
    detGb = curr.Eb .* curr.G + curr.E .* curr.Gb - 2 * curr.F .* curr.Fb;

    detG0 = orig.E .* orig.G - orig.F .^ 2;
    detG0a = orig.Ea .* orig.G + orig.E .* orig.Ga - 2 * orig.F .* orig.Fa;
    detG0b = orig.Eb .* orig.G + orig.E .* orig.Gb - 2 * orig.F .* orig.Fb;

    detC = detG ./ detG0;
    detCa = (detGa - detC .* detG0a) ./ detG0;
    detCb = (detGb - detC .* detG0b) ./ detG0;

    trC = (curr.E .* orig.G + curr.G .* orig.E - 2 * curr.F .* orig.F) ./ detG0;
    trCa = ((curr.Ea .* orig.G + curr.E .* orig.Ga - 2 * curr.Fa .* orig.F - 2 * curr.F .* orig.Fa + curr.Ga .* orig.E + curr.G .* orig.Ea) - trC .* detG0a) ./ detG0;
    trCb = ((curr.Eb .* orig.G + curr.E .* orig.Gb - 2 * curr.Fb .* orig.F - 2 * curr.F .* orig.Fb + curr.Gb .* orig.E + curr.G .* orig.Eb) - trC .* detG0b) ./ detG0;
        
%    s = sqrt(detG0);
%
%    c0 = E * (trC - 1) ./ s;
%    c0a = (E * trCa - c0 .* detG0a ./ (2 * s)) ./ s;
%    c0b = (E * trCb - c0 .* detG0b ./ (2 * s)) ./ s;
%
%    c1 = (b * (detC - 1) - E) ./ s;
%    c1a = (b * detCa - c1 .* detG0a ./ (2 * s)) ./ s;
%    c1b = (b * detCb - c1 .* detG0b ./ (2 * s)) ./ s;
    s = sqrt(detG);

    c0 = E * (trC - 1) .* s ./ detG0;
    c0a = E * trCa .* s ./ detG0 + E * (trC - 1) ./ detG0 .* (detGa ./ (2 * s) - s .* detG0a ./ detG0);
    c0b = E * trCb .* s ./ detG0 + E * (trC - 1) ./ detG0 .* (detGb ./ (2 * s) - s .* detG0b ./ detG0);

    c1 = (b * (detC - 1) - E) .* s ./ detG0;
    c1a = b * detCa .* s ./ detG0 + (b * (detC - 1) - E) ./ detG0 .* (detGa ./ (2 * s) - s .* detG0a ./ detG0);
    c1b = b * detCb .* s ./ detG0 + (b * (detC - 1) - E) ./ detG0 .* (detGb ./ (2 * s) - s .* detG0b ./ detG0);

    cE = c0 .* orig.G + c1 .* curr.G;
    cEa = c0a .* orig.G + c0 .* orig.Ga + c1a .* curr.G + c1 .* curr.Ga;
    cF = -c0 .* orig.F - c1 .* curr.F;
    cFa = -c0a .* orig.F - c0 .* orig.Fa - c1a .* curr.F - c1 .* curr.Fa;
    cFb = -c0b .* orig.F - c0 .* orig.Fb - c1b .* curr.F - c1 .* curr.Fb;
    cG = c0 .* orig.E + c1 .* curr.E;
    cGb = c0b .* orig.E + c0 .* orig.Eb + c1b .* curr.E + c1 .* curr.Eb;

    fa = [cE cE cE] .* curr.taa + [cEa cEa cEa] .* curr.ta + [cF cF cF] .* curr.tab + [cFa cFa cFa] .* curr.tb;
    fb = [cG cG cG] .* curr.tbb + [cGb cGb cGb] .* curr.tb + [cF cF cF] .* curr.tab + [cFb cFb cFb] .* curr.ta;

    f = (fa + fb) ./ [s s s];
    %f = E / 4 * ((trC - 2).^2 + 2 * (trC - 2) - 2 * (detC - 1)) + b / 4 * (detC - 1).^2
    s11 = E * (trC - 1) .* orig.G + (b * (detC - 1) - E) .* curr.G;
    s12 = E * (trC - 1) .* orig.F + (b * (detC - 1) - E) .* curr.F;
    s22 = E * (trC - 1) .* orig.E + (b * (detC - 1) - E) .* curr.E;

    data = [s11 s12 s22];
end

function [f, data] = bending_force(orig, curr, k)
    mag = k * curr.H .* (curr.H .^2 - curr.K);
    f = [mag mag mag] .* curr.n;
    data = 0;
end

function fn = property_calculator(I, dIa, dIb, dIaa, dIab, dIbb)
    function S = calculate_properties(x, y, z)
        cx = I \ x;
        cy = I \ y;
        cz = I \ z;

        ta = dIa * [cx cy cz];
        tb = dIb * [cx cy cz];

        taa = dIaa * [cx cy cz];
        tab = dIab * [cx cy cz];
        tbb = dIbb * [cx cy cz];

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

        S.cx = cx;
        S.cy = cy;
        S.cz = cz;

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
