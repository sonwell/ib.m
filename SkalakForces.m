function fn = SkalakForces(shear, bulk)
    fn = @(orig) initializer(orig, shear, bulk);
end

function fn = initializer(orig, E, b)
    oE = orig.E; oEa = orig.Ea; oEb = orig.Eb;
    oF = orig.F; oFa = orig.Fa; oFb = orig.Fb;
    oG = orig.G; oGa = orig.Ga; oGb = orig.Gb;

    detG0 = oE .* oG - oF .^ 2;
    detG0a = oEa .* oG + oE .* oGa - 2 * oF .* oFa;
    detG0b = oEb .* oG + oE .* oGb - 2 * oF .* oFb;

    s = detG0.^(-1/2);
    sa = -1/2 * detG0a .* (detG0.^(-3/2));
    sb = -1/2 * detG0b .* (detG0.^(-3/2));

    function f = forces(curr)
        cE = curr.E; cEa = curr.Ea; cEb = curr.Eb;
        cF = curr.F; cFa = curr.Fa; cFb = curr.Fb;
        cG = curr.G; cGa = curr.Ga; cGb = curr.Gb;

        detG = cE .* cG - cF .^ 2;
        detGa = cEa .* cG + cE .* cGa - 2 * cF .* cFa;
        detGb = cEb .* cG + cE .* cGb - 2 * cF .* cFb;

        detC = detG ./ detG0;
        detCa = (detGa - detC .* detG0a) ./ detG0;
        detCb = (detGb - detC .* detG0b) ./ detG0;

        trC = (cE .* oG + cG .* oE - 2 * cF .* oF) ./ detG0;
        trCa = ((cEa .* oG + cE .* oGa - 2 * cFa .* oF - 2 * cF .* oFa + ...
                 cGa .* oE + cG .* oEa) - trC .* detG0a) ./ detG0;
        trCb = ((cEb .* oG + cE .* oGb - 2 * cFb .* oF - 2 * cF .* oFb + ...
                 cGb .* oE + cG .* oEb) - trC .* detG0b) ./ detG0;

        c0 = E * (trC - 1) .* s;
        c0a = E * trCa .* s + E * (trC - 1) .* sa;
        c0b = E * trCb .* s + E * (trC - 1) .* sb;

        c1 = (b * (detC - 1) - E) .* s;
        c1a = b * detCa .* s + (b * (detC - 1) - E) .* sa;
        c1b = b * detCb .* s + (b * (detC - 1) - E) .* sb;

        coE = c0 .* oG + c1 .* cG;
        coEa = c0a .* oG + c0 .* oGa + c1a .* cG + c1 .* cGa;
        coF = (c0 .* oF + c1 .* cF);
        coFa = (c0a .* oF + c1a .* cF) + (c0 .* oFa + c1 .* cFa);
        coFb = (c0b .* oF + c0 .* oFb) + (c1b .* cF + c1 .* cFb);
        coG = c0 .* oE + c1 .* cE;
        coGb = c0b .* oE + c0 .* oEb + c1b .* cE + c1 .* cEb;

        ta = curr.ta;
        tb = curr.tb;
        taa = curr.taa;
        tab = curr.tab;
        tbb = curr.tbb;

        fa = sc(coE, taa) + sc(coEa, ta) - sc(coF, tab) - sc(coFa, tb);
        fb = sc(coG, tbb) + sc(coGb, tb) - sc(coF, tab) - sc(coFb, ta);
        f = (fa + fb) .* [s s s];
    end

    fn = @forces;
end


function v = sc(scale, u)
    v = [scale scale scale] .* u;
end
