function DM = DifferenceMatrix(dsites, ctrs)      
    [d, c] = ndgrid(dsites, ctrs);
    DM = d - c;
end