function fn = BendingForce(k)
    fn = @(obj) initializer(obj, k);
end

function fn = initializer(obj, k)
    gr = obj.geometry('reference');
    grs = gr.at_sample_sites;

    oE = grs.E; oF = grs.F; oG = grs.G;
    detG0 = oE .* oG - oF .^ 2;

    Da_d = obj.opd.d{1}; Db_d = obj.opd.d{2};
    Da_s = obj.ops.d{1}; Db_s = obj.ops.d{2};

    function f = forces()
        gc = obj.geometry('current');
        gcs = gc.at_sample_sites;
        gcd = gc.at_data_sites;

        cE = gcs.E; cF = gcs.F; cG = gcs.G;
        detG = cE .* cG - cF .^ 2;

        Hs = gcs.H; Ks = gcs.K;
        Hd = gcd.H;

        Ed = gcd.E; Fd = gcd.F; Gd = gcd.G;
        detGd = Ed .* Gd - Fd .^ 2;

        Ha_d = Da_d * Hd; Hb_d = Db_d * Hd;
        LH = (Da_s * ((Gd .* Ha_d - Fd .* Hb_d) ./ sqrt(detGd)) + ...
              Db_s * ((-Fd .* Ha_d + Ed .* Hb_d) ./ sqrt(detGd))) ./ sqrt(detG);

        detC = detG ./ detG0;
        mag = -k * (LH + Hs .* (Hs .^2 - Ks)) .* sqrt(detC);
        f = bsxfun(@times, gcs.n, mag);
    end

    fn = @forces;
end
