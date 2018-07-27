function fn = TetherForces(spring)
    fn = @(obj) initializer(obj, spring);
end

function fn = initializer(obj, spring)
    gr = obj.geometry('reference');
    grs = gr.at_sample_sites;
    ox = grs.r;
    function f = forces()
        gc = obj.geometry('current');
        gcs = gc.at_sample_sites;
        f = -spring * (gcs.r - ox);
    end

    fn = @forces;
end
