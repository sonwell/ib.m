function f = ForceCalculator(alpha, beta, gamma, delta, x0, y0, z0, youngs, bulk, bending)
    f = force_calculator(0.25, alpha, beta, gamma, delta, x0, y0, z0, youngs, bulk, bending);
end

function calculator = force_calculator(e, a, b, c, d, x0, y0, z0, youngs, bulk, bending)
    [I, Da, Db, Daa, Dab, Dbb] = RBFOperators(e, a, b, c, d);
    geometry = SurfaceGeometryCalculator(Da, Db, Daa, Dab, Dbb);
    orig = geometry(x0, y0, z0);
    
    function [X, fsk, fb, curr] = f_calculator(x, y, z)
        X = I * [x y z];
        curr = geometry(x, y, z);
        fsk = sk_force(orig, curr, youngs, bulk);
        fb = bending_force(orig, curr, bending);
    end
    calculator = @f_calculator;
end

function f = sk_force(orig, curr, E, b)
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
        
    s = sqrt(detG0);

    c0 = E * (trC - 1) ./ s;
    c0a = (E * trCa - c0 .* detG0a ./ (2 * s)) ./ s;
    c0b = (E * trCb - c0 .* detG0b ./ (2 * s)) ./ s;

    c1 = (b * (detC - 1) - E) ./ s;
    c1a = (b * detCa - c1 .* detG0a ./ (2 * s)) ./ s;
    c1b = (b * detCb - c1 .* detG0b ./ (2 * s)) ./ s;

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
end

function f = bending_force(orig, curr, k)
    detG = curr.E .* curr.G - curr.F .^ 2;
    detG0 = orig.E .* orig.G - orig.F .^ 2;
    detC = detG ./ detG0;
    mag = k * curr.H .* (curr.H .^2 - curr.K) .* sqrt(detC);
    f = [mag mag mag] .* curr.n;
end
