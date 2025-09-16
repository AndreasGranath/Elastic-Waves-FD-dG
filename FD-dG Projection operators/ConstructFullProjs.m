
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION: Assembles pairs of norm compatible good and bad projection 
%              operators from a 1D FD grid to a 1D dG grid of the same order
%              
%              accurate to a given order
%              
%
% INPUT : - x_FD: a vector containing FD nodes
%
%         - x_dG: a vector containing dG nodes
%
%         - order: the order of the dG method used 
%
% OUTPUT: - P_FD_2_dG_g / P_FD_2_dG_b: good and bad projection operators
%           from the FD grid to the dG grid.
%
%         - P_dG_2_FD_g / P_dG_2_FD_b: good and bad projection operators
%           from the dG grid to the FD grid.
%
% AUTHOR:  Andreas Granath (andreas.granath@umu.se)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [P_FD_2_dG_g,P_FD_2_dG_b,P_dG_2_FD_g,P_dG_2_FD_b]=ConstructFullProjs(x_FD,x_dG,order)
    nFD=length(x_FD);
    %ndG=length(x_dG);

    % Construct projection pairs from FD to FD-glue
      [Pf2G_b,PG2f_g,Pf2G_g,PG2f_b]=make_projection_FD(nFD,order);

   % Construct projection ops between the glues
    %[PfdG2dg,PdgG2fd]=glue_projections(AddLagrangeNodes(x_FD,order),x_dG,order);
    ndG=length(x_dG)/(order+1)+1;
    %[PdgG2fd,PfdG2dg]=NewGlueProjections(ndG,order,1);   
    [PfdG2dg,PdgG2fd]=NewGlueProjections(ndG,order,0);
    % Assemble full projection ops
    P_FD_2_dG_g=PfdG2dg*Pf2G_g;
    P_FD_2_dG_b=PfdG2dg*Pf2G_b;
    P_dG_2_FD_g=PG2f_g*PdgG2fd;
    P_dG_2_FD_b=PG2f_b*PdgG2fd;

end