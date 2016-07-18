function island = Z_island_nn(G, nn)
% 6/22/2016 Parker MacCready
%
% like Z_island but with a variable masking region
%
% DEFINE MASKING ARRAY
% island is a logical array, with:
% 1 = land (or masked out), and
% 0 = water

island = ~G.mask_rho; 
% mask out nudging regions
% ## get nn from the inputs ##
%nn = 8; % original is 8 (9 MAY make bad point go away)
%nn = 100; % sensitivity testing 6/22/2016
island(1:nn,:) = 1;
island(:,1:nn) = 1;
% and regions to north and east
island(end-1:end,:) = 1;
island(:,end-1:end) = 1;
