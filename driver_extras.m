% Driver for extra terms in the energy budget. Saves one 2D (vertically integrated)
% snapshot of extra terms, every hour. We do time averaging later.
%
% file timing and numbering:
% nn     o   nn+1
% his    o   his
%  |    nn    |
%  |  avg,dia |

% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');

% create top level output directory, if needed
odir_top = [Tdir.output,'energy_out/'];
if ~exist(odir_top,'dir'); mkdir(odir_top); end;

% get info for processing
[basename,nn_vec,dir0] = Z_runspec_raw;

% create EMPTY result directory (but don't geet rid of the flux results)
odir_base = [odir_top,basename,'/'];
if ~exist(odir_base,'dir'); mkdir(odir_base); end;
odir = [odir_top,basename,'/extras_raw/'];
if exist(odir,'dir'); rmdir(odir,'s'); end; mkdir(odir);
disp(' '); disp(['extras: Saving results to ',odir]);

%% get island mask and hypsometry

ns = num2str(nn_vec(1)); ns = ['0000',ns]; ns = ns(end-3:end);
fn = [dir0.his,'ocean_his_',ns,'.nc'];
[G,S,T] = Z_get_basic_info(fn);

% get the masking array.
% This defines the region for flattening used in APE calculation
island = Z_island(G);
G.island = island; % add island to G
save([odir,'G.mat'],'G'); % just save G in one place, to save space

% get hypsometry structure
[H] = Z_hyp(G,S,island);

%% time loop
for ii = 1:length(nn_vec);
    tic
    nn = nn_vec(ii);
    [ex2,info] = extras(dir0,nn,H,G,S);
    ncpad = ['0000',num2str(nn)]; ncpad = ncpad(end-4:end);
    save([odir,'extras_',ncpad,'.mat'],'ex2','info');
    dt = toc; disp(['file ',ncpad,' took ',num2str(round(dt)),' sec'])
end


