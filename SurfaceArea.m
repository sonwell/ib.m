function area = SurfaceArea(structure)
    id = structure.id;
    phi = structure.phi;
    npt = structure.interpolation.n;
    zero = zeros(npt, 1);
    weights = phi \ [zero; 4*pi];
    sph_sa = weights(1:npt);
    col_sums = sum(id)';
    scale = spdiags(1./col_sums, [0], npt, npt);
    esa = id * scale * sph_sa;

    transform = structure.reference;
    geometry = SurfaceGeometryCalculator(structure);
    params = structure.interpolation.parametrization;
    orig = geometry(sphere(params));
    curr = geometry(transform(params));
    da = sqrt((curr.E .* curr.G - curr.F.^2) ./ (orig.E .* orig.G - orig.F.^2));
    area = da .*  esa;
end

function x = sphere(params)
    alpha = params(:, 1);
    beta = params(:, 2);
    x = [sin(beta) .* cos(alpha), sin(beta) .* sin(alpha), cos(beta)];
end
