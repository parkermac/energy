function [H] = Z_hyp(G,S,island)
% code to get hypsometry structure H

% Note: I think the area vector (Ah) is not used by anything.

% area arrays
da = G.DX .* G.DY;
da3 = repmat(da,[1 1 S.N]); da3 = permute(da3,[3 1 2]);

% island3 is true for cells we want to mask out
island3 = repmat(island,[1 1 S.N]);
island3 = permute(island3,[3 1 2]);

% hypsometry
%
% the goal is to have functions that relates z, volume, and area.
% specifically:
% Vh(Zh) and Ah(Zh) are the volume below and area at a given z-position.
% They exist at all the depths that exist in the grid bathymetry
% (after masking using island), and with repeated depths removed.
%
% The way the indexing works is that at a given Zh the volume is
% for all cells deeper with depths deeper than that (of course the cell
% at that depth would have no volume).  In addition the area at a given
% Zh is only for cells deeper than that one, and so it does not include
% the area associated with that cell.
%
% Consistent with this indexing, the first point in the returned vectors
% has the deepest depth in the grid, zero volume, and zero area.
%
% The uppermost point is assigned an artifical z-value of 10 m, so that
% we never interpolate beyond it, and its area should be equal to the area
% of all unmasked cells.

% these become column vectors
hh = G.h(~island);
daa = da(~island);
%
[hh,ind] = sort(hh,'descend'); % deepest (greatest depth) to shallowest
daa = daa(ind); % sorted areas to go with the depths

% create three vectors, deepest to shallowest, as described above
ZZ = [-hh(2:end); 10]; % z of the next shallowest grid cell, with padding
AA = cumsum(daa); % area of all cells that are deeper than ZZ
VV = nan(size(hh)); % preallocate
for ii = 1:length(hh)
    VV(ii,1) = sum((hh(1:ii)+ZZ(ii)).*daa(1:ii));
end
% add the very bottom point
ZZ = [-hh(1); ZZ];
AA = [0; AA];
VV = [0; VV];
% VV = volume of all cells that are deeper than ZZ
%
% H.V_check = sum(daa.*hh);
% Result: V_check is exactly the same as VV(end) IF we make the "padding"
% at the top of ZZ equal to 0.  This is as expected.

%
% remove repeat values, where ZZ and VV are not changing, even though AA is
[Vh,indh] = unique(VV);
% Note: using VV seemed to get rid of a few more repeats than using ZZ

% pack results
H.Vh = Vh;
H.Zh = ZZ(indh);
H.Ah = AA(indh);
H.island = island;
H.island3 = island3;
H.da = da;
H.da3 = da3;
