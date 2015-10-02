function [D] = Z_flat(rho,fld,zr,zw,H)
% code to define a resting background state.
%
% have to run Z_hyp.m first to get hypsometry structure H
%
% The field "fld" is an extra companion to rho.  It will be processed in
% the same way as rho, but rho will be used for all sorting.
%
% The only time we use fld is for the field drho/dt, used in calculating
% the rate of change of the background APE.

%% Get some things from the hypsometry structure H
island3 = H.island3;
Vh = H.Vh;
Zh = H.Zh;
da3 = H.da3;

% make and array of grid box volumes
dz = diff(zw);
dv = dz.*da3; % array of volumes

%% now form the flattened density, and its depth profile

% form a long column vector of indices into any of the 3D matrices
% that are the same shape as rho
imat = find(ones(size(rho))==1);

% create long vectors of just the good points for all arrays we are
% going to work on
rho_vec = rho(~island3); % density
fld_vec = fld(~island3); % field
dv_vec = dv(~island3);   % volume of each grid box
imat_vec = imat(~island3); % indices

% Start by making a few long vectors that are all sorted by density:
%  density: Rf
%  field: Ff
%  grid box volume: dv_f
%  index into original arrays: Imf
% And each of these just has good points.
%
% sort on density, heaviest to lightest
[Rf,indf] = sort(rho_vec,'descend');
% then get other fields that correspond to each element in Rf
Ff = fld_vec(indf); % fld
dv_f = dv_vec(indf); % box volume
Imf = imat_vec(indf); % index into any array

% Now make a long vector which is
% volume below a given density (to within a half of a dv)
Vf = cumsum(dv_f); 

% And use the hypsometry H.Zh as a function of H.Vh to
% make z position of each flattened grid cell
Zf = interp1(Vh,Zh,Vf);
% Check: Zf(end) is within 1e-10 m of eta_flat
% disp(['Zf(end)  = ',sprintf('%.12f',Zf(end))])

% So at this point we just have long vectors of the flattened density,
% profile: Rf and Zf, (and Ff), but we need to do some interpolation
% to get this information for each grid box.

% Make the flattened density at each in situ z position.
% i.e. make Rf at zr.
% Eventually we will use this
% First add some padding before interpolating
% NOTE: changing the -500 to -10 made little difference
Zf_ext = [Zf(1) - 500; Zf; Zf(end) + 10];
Rf_ext = [Rf(1); Rf; Rf(end)];
% Rf_at_z is the flattened density at every good zr position.
% Note that zr(Imf) is a vector of the original z position of
% every good grid box, sorted by density.
Rf_at_z = interp1(Zf_ext,Rf_ext,zr(Imf));
Ff_ext = [Ff(1); Ff; Ff(end)];
Ff_at_z = interp1(Zf_ext,Ff_ext,zr(Imf));
% and these are all vectors, with length equal to the number of good
% elements in the original fields

% Now make some vertical integrals
dZf = diff(Zf);
dZf_ext = diff(Zf_ext);
% because of the way these are used it doesn't matter that they
% start from the bottom of the "ext" depth

% NOTE: this integral seems suspicious
% START HERE
intRf = cumsum(Rf.*dZf_ext(1:end-1));
intRf_ext = [0; intRf; intRf(end)];
    
intRf_at_z = interp1(Zf_ext,intRf_ext,zr(Imf));
% compare to the line from a later section:
% intRf_at_zf = interp1(Zf_ext,intRf_ext,Zf_mat(Imf));

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

% initializing
R_mat = rho;
F_mat = fld;
Z_mat = zr;
DV_mat = dv;
% masking
R_mat(island3) = NaN;
F_mat(island3) = NaN;
Z_mat(island3) = NaN;
DV_mat(island3) = NaN;
% filling flattened fields
nmat = nan(size(R_mat));
Rf_mat = nmat;
Rf_mat(Imf) = Rf_at_z;
Ff_mat = nmat;
Ff_mat(Imf) = Ff_at_z;
Zf_mat = nmat;
Zf_mat(Imf) = Zf;

%% save results to a structure, with shorter names

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
