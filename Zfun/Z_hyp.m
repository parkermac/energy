function [H] = Z_hyp(G,S,island)
% 6/12/2014  Parker MacCready
% code to get hypsometry structure H

% area and volume arrays
da = G.DX .* G.DY;
da3 = repmat(da,[1 1 S.N]); da3 = permute(da3,[3 1 2]);

% island3 is true for cells we want to mask out
island3 = repmat(island,[1 1 S.N]);
island3 = permute(island3,[3 1 2]);

% hypsometry
%
% the goal is to have a function that relates z, volume, and area
% in practice we will need:
% 1. z = f(volume)
% 2. area = ff(z) = d(volume)/dz
%
% these become column vectors
hh = G.h(~island);
daa = da(~island);
%
[hh,ind] = sort(hh,'descend'); % deepest to shallowest
daa = daa(ind);
% three vectors, deepest to shallowest, relating z, volume below z,
% and area at z
ZZ = [-hh(2:end); 10]; % z of the next shallowest grid cell, with padding
AA = cumsum(daa); % area of all cells that are deeper than ZZ
%
for ii = 1:length(hh)
    VV(ii,1) = sum((hh(1:ii)+ZZ(ii)).*daa(1:ii));
end
% add the very bottom point
ZZ = [-hh(1); ZZ];
AA = [0; AA];
VV = [0; VV];
% VV = volume below ZZ of all cells that are deeper than ZZ
%
% V_check = sum(daa.*(hh + eta_flat));
% Result: V_check is exactly the same as VV(end)
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
