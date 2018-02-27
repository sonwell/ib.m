function X = RunTest(domain, shape)
    load nodesets/sphere.mat;
    X = RBFFluidTest(domain, sphere.data_sites, sphere.sample_sites, 1, 8.9e-4, shape);
end
