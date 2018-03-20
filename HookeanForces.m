function fn = HookeanForces(spring)
    fn = @(orig) initializer(orig, spring);
end

function fn = initializer(orig, E)
    ox = orig.x;
    function f = forces(curr)
        f = -E * (curr.x - ox);
    end

    fn = @forces;
end
