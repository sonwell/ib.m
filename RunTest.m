function X = RunTest(domain)
    cell = RedBloodCell(100, 400, SkalakForces(2.5e-3, 2.5e-1));
    vessel = BloodVessel(400, 1600, HookeanForces(1e6));
    X = RBFFluidTest(domain, 1, 8.9e-3, cell, vessel);
end
