function [D] = Z_flat_fast(rho,zr,zw,H)
% 7/19/2015  Parker MacCready
% code to define a resting background state.
%
% have to run Z_hyp.m first to get hypsometry structure H
%
% Omits the second field used by Z_flat (only about a 15% speedup).

island3 = H.island3;
Vh = H.Vh;
Zh = H.Zh;
da3 = H.da3;

dz = diff(zw);
dv = dz.*da3; % array of volumes

%% now form the flattened density, and Winters's z_star

% form an index matrix
imat = find(ones(size(rho))==1);

rho_vec = rho(~island3);
dv_vec = dv(~island3);
imat_vec = imat(~island3);

[Rf,indf] = sort(rho_vec,'descend');

dv_f = dv_vec(indf); % same as Vh(end) to 3 in 1e14.

Imf = imat_vec(indf); % the index into rho where each Zf goes. Nice!

Vf = cumsum(dv_f); % volume below a given density (given by Rf)

Zf = interp1(Vh,Zh,Vf);
% Check Zf(end) is deeper than eta_flat by 2e-12 m! Good

% make Rf at z
% add some padding before interpolating
% NOTE: changing the -500 to -10 made little difference
Zf_ext = [Zf(1) - 500; Zf; Zf(end) + 10];
Rf_ext = [Rf(1); Rf; Rf(end)];
Rf_at_z = interp1(Zf_ext,Rf_ext,zr(Imf));

% make the integral required for more-perfect APE
% Scotti & White, but incompressible

dZf_ext = diff(Zf_ext);
intRf = cumsum(Rf.*dZf_ext(1:end-1));
intRf_ext = [0; intRf; intRf(end)];

intRf_at_z = interp1(Zf_ext,intRf_ext,zr(Imf));
% compare to the line from a later section:
% intRf_at_zf = interp1(Zf_ext,intRf_ext,Zf_mat(Imf));

%% so now we have the means to form 3D arrays of all the needed fields
%
% we already have rho and z, and now using the index Imf we add the
% fields Zf (the z_flat for this rho), and Rf (the value of rho_flat
% at this z).

R_mat = rho;
Z_mat = zr;
DV_mat = dv;
R_mat(island3) = NaN;
Z_mat(island3) = NaN;
DV_mat(island3) = NaN;
nmat = nan(size(R_mat));
Rf_mat = nmat;
Rf_mat(Imf) = Rf_at_z;
Zf_mat = nmat;
Zf_mat(Imf) = Zf;

%% save results to a structure
D.R = R_mat;
D.Z = Z_mat;
D.DV = DV_mat;
D.Rf = Rf_mat;
D.Zf = Zf_mat;

intRf_mat = nmat;
intRf_mat(Imf) = intRf_at_z;
intRf_at_zf = interp1(Zf_ext,intRf_ext,Zf_mat(Imf));
intRfzf_mat = nmat;
intRfzf_mat(Imf) = intRf_at_zf;
D.intRf = intRf_mat;
D.intRfzf = intRfzf_mat;

