function fn = SkalakForces(shear, bulk)
    fn = @(obj) initializer(obj, shear, bulk);
end

function fn = initializer(obj, shear, bulk)
    gr = obj.geometry('reference');
    grs = gr.at_sample_sites;

    oE = grs.E; oEa = 2 * dot(grs.dr{1}, grs.ddr{1, 1}, 2); oEb = 2 * dot(grs.dr{1}, grs.ddr{1, 2}, 2);
    oG = grs.G; oGa = 2 * dot(grs.dr{2}, grs.ddr{1, 2}, 2); oGb = 2 * dot(grs.dr{2}, grs.ddr{2, 2}, 2);
    oF = grs.F; oFa = dot(grs.dr{1}, grs.ddr{1, 2}, 2) + dot(grs.dr{2}, grs.ddr{1, 1}, 2);
                oFb = dot(grs.dr{1}, grs.ddr{2, 2}, 2) + dot(grs.dr{2}, grs.ddr{1, 2}, 2);

    detG0 = oE .* oG - oF .^ 2;
    detG0a = oEa .* oG + oE .* oGa - 2 * oF .* oFa;
    detG0b = oEb .* oG + oE .* oGb - 2 * oF .* oFb;

    s = detG0.^(-1/2);
    sa = -1/2 * detG0a .* (detG0.^(-3/2));
    sb = -1/2 * detG0b .* (detG0.^(-3/2));

    function f = forces()
        gc = obj.geometry('current');
        gcs = gc.at_sample_sites;

        cE = gcs.E; cEa = 2 * dot(gcs.dr{1}, gcs.ddr{1, 1}, 2); cEb = 2 * dot(gcs.dr{1}, gcs.ddr{1, 2}, 2);
        cG = gcs.G; cGa = 2 * dot(gcs.dr{2}, gcs.ddr{1, 2}, 2); cGb = 2 * dot(gcs.dr{2}, gcs.ddr{2, 2}, 2);
        cF = gcs.F; cFa = dot(gcs.dr{1}, gcs.ddr{1, 2}, 2) + dot(gcs.dr{2}, gcs.ddr{1, 1}, 2);
                    cFb = dot(gcs.dr{1}, gcs.ddr{2, 2}, 2) + dot(gcs.dr{2}, gcs.ddr{1, 2}, 2);

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

        c0 = shear * (trC - 1) .* s;
        c0a = shear * trCa .* s + shear * (trC - 1) .* sa;
        c0b = shear * trCb .* s + shear * (trC - 1) .* sb;

        c1 = (bulk * (detC - 1) - shear) .* s;
        c1a = bulk * detCa .* s + (bulk * (detC - 1) - shear) .* sa;
        c1b = bulk * detCb .* s + (bulk * (detC - 1) - shear) .* sb;

        coE = c0 .* oG + c1 .* cG;
        coEa = c0a .* oG + c0 .* oGa + c1a .* cG + c1 .* cGa;
        coF = (c0 .* oF + c1 .* cF);
        coFa = (c0a .* oF + c1a .* cF) + (c0 .* oFa + c1 .* cFa);
        coFb = (c0b .* oF + c0 .* oFb) + (c1b .* cF + c1 .* cFb);
        coG = c0 .* oE + c1 .* cE;
        coGb = c0b .* oE + c0 .* oEb + c1b .* cE + c1 .* cEb;

        ta = gcs.dr{1};
        tb = gcs.dr{2};
        taa = gcs.ddr{1, 1};
        tab = gcs.ddr{1, 2};
        tbb = gcs.ddr{2, 2};

        sc = @(a, b) bsxfun(@times, b, a);

        fa = sc(coE, taa) + sc(coEa, ta) - sc(coF, tab) - sc(coFa, tb);
        fb = sc(coG, tbb) + sc(coGb, tb) - sc(coF, tab) - sc(coFb, ta);
        f = sc(s, fa + fb);
    end

    fn = @forces;
end
