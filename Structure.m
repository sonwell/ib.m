function structure = Structure(params, reference, varargin)
    interpolation.parametrization = params.data_sites;
    interpolation.n = size(params.data_sites, 1);

    evaluation.parametrization = params.sample_sites;
    evaluation.n = size(params.sample_sites, 1);

    [id, da, db, daa, dab, dbb, phi] = RBFOperators(params.data_sites, params.sample_sites);
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

    structure.force = CombineForces(varargin{:});
end
