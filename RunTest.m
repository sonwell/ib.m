function X = RunTest(domain)
    cell = RedBloodCell(625, 3025, SkalakForces(2.5e-3, 2.5e-1));
    vessel = BloodVessel(1600, 12800, HookeanForces(2.45e-3)); %1e8));
    X = RBFFluidTest(domain, 1, 8.9e-3, cell, vessel);
end
