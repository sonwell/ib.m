function structure = Structure(eps, alpha, beta, gamma, delta, reference)
    [id, da, db, daa, dab, dbb, phi] = RBFOperators(eps, alpha, beta, alpha, beta);
    interpolation.lambda = beta;
    interpolation.theta = alpha;
    interpolation.n = size(alpha, 1);

    [id, da, db, daa, dab, dbb, phi] = RBFOperators(eps, alpha, beta, gamma, delta);
    evaluation.lambda = delta;
    evaluation.theta = gamma;
    evaluation.n = size(gamma, 1);

    structure.eps = eps;
    structure.interpolation = interpolation;
    structure.evaluation = evaluation;
    structure.reference = reference;
    structure.id = id;
    structure.da = da;
    structure.db = db;
    structure.daa = daa;
    structure.dab = dab;
    structure.dbb = dbb;
    structure.phi = phi;
end
