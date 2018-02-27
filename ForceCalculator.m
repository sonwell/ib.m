function f = ForceCalculator(structure, youngs, bulk, bending)
    f = force_calculator(structure, youngs, bulk, bending);
end

function calculator = force_calculator(structure, youngs, bulk, bending)
    transform = structure.reference;
    alpha = structure.interpolation.theta;
    beta = structure.interpolation.lambda;
    id = structure.id;
    
    X0 = transform(alpha, beta);
    geometry = SurfaceGeometryCalculator(structure);
    orig = geometry(X0);

    function [f, x] = f_calculator(X)
        curr = geometry(X);
        fsk = sk_force(orig, curr, youngs, bulk);
        fb = bending_force(orig, curr, bending);
        f = fsk + fb;
        x = id * X;
    end
    calculator = @f_calculator;
end

function f = sk_force(orig, curr, E, b)
    oE = orig.E; cE = curr.E;
    oF = orig.F; cF = curr.F;
    oG = orig.G; cG = curr.G;

    oEa = orig.Ea; cEa = curr.Ea;
    oFa = orig.Fa; cFa = curr.Fa;
    oGa = orig.Ga; cGa = curr.Ga;

    oEb = orig.Eb; cEb = curr.Eb;
    oFb = orig.Fb; cFb = curr.Fb;
    oGb = orig.Gb; cGb = curr.Gb;

    detG = cE .* cG - cF .^ 2;
    detGa = cEa .* cG + cE .* cGa - 2 * cF .* cFa;
    detGb = cEb .* cG + cE .* cGb - 2 * cF .* cFb;

    detG0 = oE .* oG - oF .^ 2;
    detG0a = oEa .* oG + oE .* oGa - 2 * oF .* oFa;
    detG0b = oEb .* oG + oE .* oGb - 2 * oF .* oFb;

    detC = detG ./ detG0;
    detCa = (detGa - detC .* detG0a) ./ detG0;
    detCb = (detGb - detC .* detG0b) ./ detG0;

    trC = (cE .* oG + cG .* oE - 2 * cF .* oF) ./ detG0;
    trCa = ((cEa .* oG + cE .* oGa - 2 * cFa .* oF - 2 * cF .* oFa + cGa .* oE + cG .* oEa) - trC .* detG0a) ./ detG0;
    trCb = ((cEb .* oG + cE .* oGb - 2 * cFb .* oF - 2 * cF .* oFb + cGb .* oE + cG .* oEb) - trC .* detG0b) ./ detG0;

    s = detG0.^(-1/2);
    sa = -1/2 * detG0a .* (detG0.^(-3/2));
    sb = -1/2 * detG0b .* (detG0.^(-3/2));

    c0 = E * (trC - 1) .* s;
    c0a = E * trCa .* s + E * (trC - 1) .* sa;
    c0b = E * trCb .* s + E * (trC - 1) .* sb;

    c1 = (b * (detC - 1) - E) .* s;
    c1a = b * detCa .* s + (b * (detC - 1) - E) .* sa;
    c1b = b * detCb .* s + (b * (detC - 1) - E) .* sb;

    coE = c0 .* oG + c1 .* cG;
    coEa = c0a .* oG + c0 .* oGa + c1a .* cG + c1 .* cGa;
    coF = c0 .* oF + c1 .* cF;
    coFa = (c0a .* oF + c1a .* cF) + (c0 .* oFa + c1 .* cFa);
    coFb = (c0b .* oF + c0 .* oFb) + (c1b .* cF + c1 .* cFb);
    coG = c0 .* oE + c1 .* cE;
    coGb = c0b .* oE + c0 .* oEb + c1b .* cE + c1 .* cEb;

    fa = [coE coE coE] .* curr.taa + [coEa coEa coEa] .* curr.ta - [coF coF coF] .* curr.tab - [coFa coFa coFa] .* curr.tb;
    fb = [coG coG coG] .* curr.tbb + [coGb coGb coGb] .* curr.tb - [coF coF coF] .* curr.tab - [coFb coFb coFb] .* curr.ta;

    f = (fa + fb) .* [s s s];
end

function f = bending_force(orig, curr, k)
    detG = curr.E .* curr.G - curr.F .^ 2;
    detG0 = orig.E .* orig.G - orig.F .^ 2;
    detC = detG ./ detG0;
    mag = k * curr.H .* (curr.H .^2 - curr.K) .* sqrt(detC);
    f = [mag mag mag] .* curr.n;
end
