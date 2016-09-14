load nodesets/sphere4.mat;
try
    RBFFluidTest(20, sphere.data_sites, sphere.sample_sites, 4.096e-9, 1.6e-5);
catch E
    E
    quit()
end
