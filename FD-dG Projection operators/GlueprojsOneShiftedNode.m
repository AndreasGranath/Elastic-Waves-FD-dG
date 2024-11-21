function [PfdG2dg,PdgG2fd]=GlueprojsOneShiftedNode(xFD,xdG)
    order=3;
    g=uniquetol([xFD xdG], 10*eps);
    [PFD2g,Pg2FD]=FineToCoarseProj(xFD,g,order);
    [PdG2g,Pg2dG]=FineToCoarseProj(xdG,g,order);

    PfdG2dg=Pg2dG*PFD2g;
    PdgG2fd=Pg2FD*PdG2g;
end