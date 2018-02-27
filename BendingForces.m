function fn = BendingForce(k)
    fn = @(orig) initializer(orig, k);
end

function fn = initializer(orig, k)
    detG0 = orig.E .* orig.G - orig.F .^ 2;

    % TODO Add Laplace-Beltrami term
    function f = forces(curr)
        detG = curr.E .* curr.G - curr.F .^ 2;
        detC = detG ./ detG0;
        mag = k * curr.H .* (curr.H .^2 - curr.K) .* sqrt(detC);
        f = [mag mag mag] .* curr.n;
    end

    fn = @forces;
end
