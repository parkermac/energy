function [D] = Z_flat(eta,rho,fld,zr,zw,H)
% 5/5/2015  Parker MacCready
% code to define a resting background state.
%
% have to run Z_hyp.m first to get hypsometry structure H
%
% The field "fld" is an extra companion to rho.  It will be processed in
% the same way as rho, but rho will be used for all sorting.

%island = H.island;
island3 = H.island3;
Vh = H.Vh;
Zh = H.Zh;
%Ah = H.Ah;
%da = H.da;
da3 = H.da3;

dz = diff(zw);
dv = dz.*da3; % array of volumes

%% now form the flattened density, and Winters's z_star

% flattened surface height (not used?)
%eta_flat = sum(eta(~island).*da(~island))/sum(da(~island));

% form an index matrix
imat = find(ones(size(rho))==1);

rho_vec = rho(~island3);
fld_vec = fld(~island3);
dv_vec = dv(~island3);
imat_vec = imat(~island3);

[Rf,indf] = sort(rho_vec,'descend');
Ff = fld_vec(indf); % fld, sorted by rho

dv_f = dv_vec(indf); % same as Vh(end) to 3 in 1e14.

Imf = imat_vec(indf); % the index into rho where each Zf goes. Nice!

Vf = cumsum(dv_f); % volume below a given density (given by Rf)

%Af = interp1(Vh,Ah,Vf); % area at top of that volume
% Check: Af and Vf form essentially the same curve as Ah and Vh.

Zf = interp1(Vh,Zh,Vf);
% Check Zf(end) is deeper than eta_flat by 2e-12 m! Good

% make Rf at z
% add some padding before interpolating
Zf_ext = [Zf(1) - 500; Zf; Zf(end) + 10];
Rf_ext = [Rf(1); Rf; Rf(end)];
Rf_at_z = interp1(Zf_ext,Rf_ext,zr(Imf));
Ff_ext = [Ff(1); Ff; Ff(end)];
Ff_at_z = interp1(Zf_ext,Ff_ext,zr(Imf));

% make the integral required for more-perfect APE
% Scotti & White, but incompressible
dZf = diff(Zf);
intRf = cumsum(Rf(2:end).*dZf);
intRf = [0; intRf]; % a vector the same length as Zf
intRf_ext = [intRf(1); intRf; intRf(end)];
intRf_at_z = interp1(Zf_ext,intRf_ext,zr(Imf));
% and now for fld
intFf = cumsum(Ff(2:end).*dZf);
intFf = [0; intFf]; % a vector the same length as Zf
intFf_ext = [intFf(1); intFf; intFf(end)];
intFf_at_z = interp1(Zf_ext,intFf_ext,zr(Imf));


%% so now we have the means to form 3D arrays of all the needed fields
%
% we already have rho and z, and now using the index Imf we add the
% fields Zf (the z_flat for this rho), and Rf (the value of rho_flat
% at this z).

R_mat = rho;
F_mat = fld;
Z_mat = zr;
DV_mat = dv;
R_mat(island3) = NaN;
F_mat(island3) = NaN;
Z_mat(island3) = NaN;
DV_mat(island3) = NaN;
nmat = nan(size(R_mat));
Rf_mat = nmat;
Rf_mat(Imf) = Rf_at_z;
Ff_mat = nmat;
Ff_mat(Imf) = Ff_at_z;
Zf_mat = nmat;
Zf_mat(Imf) = Zf;

%% save results to a structure
D.R = R_mat;
D.F = F_mat;
D.Z = Z_mat;
D.DV = DV_mat;
D.Rf = Rf_mat;
D.Ff = Ff_mat;
D.Zf = Zf_mat;

intRf_mat = nmat;
intRf_mat(Imf) = intRf_at_z;
intRf_at_zf = interp1(Zf_ext,intRf_ext,Zf_mat(Imf));
intRfzf_mat = nmat;
intRfzf_mat(Imf) = intRf_at_zf;
D.intRf = intRf_mat;
D.intRfzf = intRfzf_mat;

intFf_mat = nmat;
intFf_mat(Imf) = intFf_at_z;
intFf_at_zf = interp1(Zf_ext,intFf_ext,Zf_mat(Imf));
intFfzf_mat = nmat;
intFfzf_mat(Imf) = intFf_at_zf;
D.intFf = intFf_mat;
D.intFfzf = intFfzf_mat;
