% Code to test the effect of compressibility in the calculation of APE.

% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');

% choose the file to work on
basename = 'D2005';
nsave = 1848; %
dir1 = '/Users/PM5/Documents/roms/output/';
dir0 = [dir1,basename,'_avg/'];
ns = num2str(nsave); ns = ['0000',ns]; ns = ns(end-3:end);
fn = [dir0,'ocean_avg_',ns,'.nc'];

% lat, lon indices for profile
nr = 30; nc = 30;

% depth range for APE calculation
deltaz = 500;

% constants
g = 9.81;

%% get info and 3D fields

[G,S,T] = Z_get_basic_info(fn);

[z_avg,z_w_avg] = Z_s2z(G.h,0*G.h,S);

salt = nc_varget(fn,'salt');
temp = nc_varget(fn,'temp');

[den,den1,alpha,beta] = Z_roms_eos_vec(salt,temp,z_avg);

%% get values at a single profile

SS = squeeze(salt(:,nr,nc));
TT = squeeze(temp(:,nr,nc));

DD = squeeze(den(:,nr,nc));
DD1 = squeeze(den1(:,nr,nc));
Z = squeeze(z_avg(:,nr,nc));

% interpolate to get more points
z = linspace(Z(1),Z(end),1000);
ss = interp1(Z,SS,z);
tt = interp1(Z,TT,z);
dd = interp1(Z,DD,z);
dd1 = interp1(Z,DD1,z);
dz = diff(z); dz = [dz(1), dz];

%% get values over a single range, for plotting
zdeep = -deltaz - 6; zshallow = zdeep + deltaz;

n0 = find(z>=zdeep); n0 = n0(1);
n1 = find(z>=zshallow); n1 = n1(1);

sr = ss(n0:n1); tr = tt(n0:n1);
zr = z(n0:n1); dzr = dz(n0:n1);
ddr = dd(n0:n1); dd1r = dd1(n0:n1);

ovec = ones(size(sr));

% get densities for the parcel at the bottom of the range
[dpr,dp1r,alphar,betar] = Z_roms_eos_vec(sr(1)*ovec,tr(1)*ovec,zr);

Dd = dpr-ddr;
Dd1 = dp1r-dd1r;

% calculate APE per unit area (J m-2)
dzr_cut = dzr; dzr_cut(1) = dzr(1)/2; dzr_cut(end) = dzr(end)/2;
ape = sum(g*Dd.*dzr_cut);
ape1 = sum(g*Dd1.*dzr_cut);

%% plotting, part 1

close all
lw = 3;

figure
Z_fig(16);
set(gcf,'position',[10 10 2000 1000])
set(gcf,'PaperPositionMode','auto');


subplot(131)
plot(dd,z,'-r',dd1,z,'-b','linewidth',lw);
hold on
plot(dpr,zr,'-r',dp1r,zr,'-b','linewidth',lw);
xlabel('Density Anomaly (kg m^{-3})');
ylabel('Z (m)');
title('Full Water Column')

subplot(132)
plot(Dd,zr,'-r',Dd1,zr,'-b','linewidth',lw)
xlabel('\Delta Density Anomaly (kg m^{-3})');
ylabel('Z (m)');
title([num2str(round(zr(1))),' m to ',num2str(round(zr(end))),' m']);

aa = axis;
dx = aa(2) - aa(1);
dy = aa(4) - aa(3);
text(aa(2) - .05*dx, aa(3) + .05*dy, ...
    ['Compressible APE_{V} = ',num2str(round(ape)),' J m^{-3}'], ...
    'horizontalalignment','r','color','r');
text(aa(2) - .05*dx, aa(3) + .1*dy, ...
    ['Approximate APE_{V} = ',num2str(round(ape1)),' J m^{-3}'], ...
    'horizontalalignment','r','color','b');

%% get values over a series of ranges

dzd = 50;
zdeep_vec = [-2950:dzd:-deltaz-10];
zshallow_vec = zdeep_vec + deltaz;
ape_vec = nan*zdeep_vec; % compressible
ape1_vec = ape_vec; % incompressible

for nn = 1:length(zdeep_vec)

    zdeep = zdeep_vec(nn);
    zshallow = zshallow_vec(nn);
    
    n0 = find(z>=zdeep); n0 = n0(1);
    n1 = find(z>=zshallow); n1 = n1(1);

    sr = ss(n0:n1); tr = tt(n0:n1);
    zr = z(n0:n1); dzr = dz(n0:n1);
    ddr = dd(n0:n1); dd1r = dd1(n0:n1);

    ovec = ones(size(sr));

    % get densities for the parcel at the bottom of the range
    [dpr,dp1r,alphar,betar] = Z_roms_eos_vec(sr(1)*ovec,tr(1)*ovec,zr);

    Dd = dpr-ddr;
    Dd1 = dp1r-dd1r;

    % calculate APE per unit area (J m-2)
    dzr_cut = dzr; dzr_cut(1) = dzr(1)/2; dzr_cut(end) = dzr(end)/2;
    ape_vec(nn) = sum(g*Dd.*dzr_cut);
    ape1_vec(nn) = sum(g*Dd1.*dzr_cut);
    
end

%% plotting, part 2

subplot(133)
plot(ape_vec,zshallow_vec,'-r',ape1_vec,zshallow_vec,'-b','linewidth',lw)
xlabel(['APE_{V} (J m^{-3}) for ',num2str(deltaz),' m Displacement']);
ylabel('Starting Z (m)');
ape_a = sum(dzd*ape_vec);
ape1_a = sum(dzd*ape1_vec);

aa = axis;
dx = aa(2) - aa(1);
dy = aa(4) - aa(3);
text(aa(2) - .05*dx, aa(3) + .05*dy, ...
    ['APE_{A} = ',num2str(round(ape_a/1e3)),' kJ m^{-2}'], ...
    'horizontalalignment','r','color','r');
text(aa(2) - .05*dx, aa(3) + .1*dy, ...
    ['APE_{A} = ',num2str(round(ape1_a/1e3)),' kJ m^{-2}'], ...
    'horizontalalignment','r','color','b');
text(aa(2) - .05*dx, aa(3) + .15*dy, ...
    ['Error = ',num2str(100*(ape1_a/ape_a)-100),'%'], ...
    'horizontalalignment','r','color','k');








