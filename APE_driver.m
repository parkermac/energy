function APE_driver(calc_name)
%% Driver for calculating APE of arbitrary runs (for SciDAC).
%
% Calculates the exact APE (but neglecting compressibility).

% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');

% get info for processing
switch calc_name
    case 'mac0' % single time, for testing
        basename = 'D2005';
        nn_vec = 1836; %
        dir00 = '/Users/PM5/Documents/roms/output/';
    case 'mac1' % a longer series, suitable for averaging
        basename = 'D2005';
        nn_vec = 1800:1928; % avg & dia - have 1800:1928 (his to 1929)
        dir00 = '/Users/PM5/Documents/roms/output/';
    case 'fjo' % fjord, for the full year
        basename = 'Cdia2005';
        nn_vec = 1:8760; % avg & dia - have 1:8760 (his to 8761)
        dir00 = '/data1/parker/roms/output/';
    case 's45_2005' 
        basename = 'S45_OcAcT2005R2005_2005';
        nn_vec = 1:366; 
        dir00 = '/pmraid4/parker/roms/output/';
    case 's45_2099' 
        basename = 'S45_OcAcT2005R2005_2099';
        nn_vec = 1:366; 
        dir00 = '/pmraid4/parker/roms/output/';
    case 's85_2005' 
        basename = 'S85_OcAcT2005R2005_2005';
        nn_vec = 1:366; 
        dir00 = '/pmraid4/parker/roms/output/';
    case 's85_2099' 
        basename = 'S85_OcAcT2005R2005_2099';
        nn_vec = 1:366; 
        dir00 = '/pmraid4/parker/roms/output/';        
end
dir0 = [dir00,basename,'_his/'];

% create top level output directory, if needed
odir_top = [Tdir.output,'APE_out/'];
if ~exist(odir_top,'dir'); mkdir(odir_top); end;

% create EMPTY result directory
odir = [odir_top,basename,'/'];
if exist(odir,'dir'); rmdir(odir,'s'); end;
mkdir(odir);
disp(' '); disp(['APE: Saving results to ',odir]);

%% Initializing: get island mask and hypsometry
ns = num2str(nn_vec(1)); ns = ['0000',ns]; ns = ns(end-3:end);
fn = [dir0,'ocean_his_',ns,'.nc'];
[G,S,T] = Z_get_basic_info(fn);
DA.initial = G.DX .* G.DY;

% Define the region for flattening used in APE calculation
island_temp = Z_island(G);

% get hypsometry structure
[H] = Z_hyp(G,S,island_temp);

% initialize output arrays
island_tag_list = {'full','salish','shelf'};

% make a more restricted version of the mask
sal = load('Zfun/mask_salish.mat'); % has vectors x and y
insalish = inpolygon(G.lon_rho,G.lat_rho,sal.x,sal.y);

for ii = 1:length(island_tag_list)
    island_tag = island_tag_list{ii};
    switch island_tag
        case 'full'
            island_temp = Z_island(G);
        case 'salish'
            island_temp = Z_island(G);
            island_temp(~insalish) = 1;
        case 'shelf'
            island_temp = Z_island(G);
            island_temp(G.h>400) = 1;
            island_temp(insalish) = 1;
    end

    island_temp = logical(island_temp);
    DA.(island_tag) = DA.initial;
    DA.(island_tag)(island_temp) = NaN;
    
    ape_ser.(island_tag) = nan(length(nn_vec),1);
    
end
ape_ser.td = nan(length(nn_vec),1);


%% time loop
for nn = 1:length(nn_vec);
    tic
    nnn = nn_vec(nn);
    ns = num2str(nnn); ns = ['0000',ns]; ns = ns(end-3:end);
    fn = [dir0,'ocean_his_',ns,'.nc'];
    
    [T] = Z_get_time_info(fn);
    ape_ser.td(nn) = T.time_datenum;
    
    [G, apea, apea_up, apea_down, zz, rr] = Z_ape(fn);
    
    for ii = 1:length(island_tag_list)
        island_tag = island_tag_list{ii};
        ape_ser.(island_tag)(nn) = squeeze(nansum(apea(:).*DA.(island_tag)(:)));
    end
    dt = toc;
    disp(['Done with save number ',num2str(nnn),' (',num2str(dt),' sec)']);
   
end

save([odir,'APE.mat'],'ape_ser');


