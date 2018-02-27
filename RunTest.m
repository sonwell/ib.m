function X = RunTest(domain)
    data = load('nodesets/sphere.mat', 'sphere');
    structure = Structure(data.sphere, @RBCShape, SkalakForces(2.5e-3, 2.5e-1));
    X = RBFFluidTest(domain, 1, 8.9e-4, structure);
end
