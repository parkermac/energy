function island = Z_island(G)
% 6/12/2014 Parker MacCready
%
% DEFINE MASKING ARRAY
% island is a logical array, with:
% 1 = land (or masked out), and
% 0 = water

island = ~G.mask_rho; 
% mask out nudging regions
nn = 8; % original is 8 (9 MAY make bad point go away)
island(1:nn,:) = 1;
island(:,1:nn) = 1;
% and regions to north and east
island(end-1:end,:) = 1;
island(:,end-1:end) = 1;
