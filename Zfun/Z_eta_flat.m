function eta_flat = Z_eta_flat(eta, da, island)
% gives the (scalar) value of the flattened surface height
%
% all three inputs are 2D arrays

eta_flat = sum(eta(~island).*da(~island))/sum(da(~island));