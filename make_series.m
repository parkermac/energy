% make_series.m  6/18/2014  Parker MacCready
%

clear
addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');

if 0
    pth = '/Users/PM5/Documents/tools_output/energy_out/Cdia2005/flux_raw/';
    nn_vec = 1:8760;
else
    pth = '/Users/PM5/Documents/tools_output/energy_out/Cdia2005/flux_lp71/';
    nn_vec = 36:24:8760-36;
end

island_tag = 'full'; % salish, shelf, abyss, full

load([pth,'G.mat']); % has island

td_vec = NaN * nn_vec;

tt = 0; % a counter
for nn = nn_vec
    
    tt = tt+1;
    
    ncpad = ['0000',num2str(nn)]; ncpad = ncpad(end-4:end);
    infile = ['flux_',ncpad,'.mat'];
    if(mod(nn,10)==0); disp(['Working on ',infile]); end;
    load([pth,infile]);
    td_vec(tt) = info.td_avg;
    [p2,k2,e2,i2,sw,info] = Z_make_derived(p2,k2,sw,info);
    
    if tt==1
        
        DA = G.DX .* G.DY;

        % get rid of bad river points (instead of G.island)
        island = isnan(k2.accel_check);
        
        % make a more restricted version of the mask
        sal = load('Zfun/mask_salish.mat'); % has vectors x and y
        insalish = inpolygon(G.lon_rho,G.lat_rho,sal.x,sal.y);
        switch island_tag
            case 'full'
                % no other changes
            case 'shelf_and_salish'
                island(G.h>400) = 1;
            case 'salish'
                island(~insalish) = 1;
            case 'shelf'
                island(G.h>400) = 1;
                island(insalish) = 1;
            case 'abyss'
                island(G.h<=400) = 1;
            case 'shelf_and_abyss'
                island(insalish) = 1;
        end
        
        island = logical(island);
        DA(island) = NaN;
        
    end
    
    % do area integrals to get volume integrals (Watts)
    % (or Joules for ape_avg)
    
    fnm = fieldnames(p2);
    for ii = 1:length(fnm)
        varname = fnm{ii};
        P2.(varname)(tt) = squeeze(nansum(p2.(varname)(:).*DA(:)));
    end
    fnm = fieldnames(k2);
    for ii = 1:length(fnm)
        varname = fnm{ii};
        K2.(varname)(tt) = squeeze(nansum(k2.(varname)(:).*DA(:)));
    end
    fnm = fieldnames(e2);
    for ii = 1:length(fnm)
        varname = fnm{ii};
        E2.(varname)(tt) = squeeze(nansum(e2.(varname)(:).*DA(:)));
    end
    fnm = fieldnames(i2);
    for ii = 1:length(fnm)
        varname = fnm{ii};
        I2.(varname)(tt) = squeeze(nansum(i2.(varname)(:).*DA(:)));
    end
    fnm = fieldnames(sw);
    for ii = 1:length(fnm)
        varname = fnm{ii};
        SW.(varname)(tt) = squeeze(nansum(sw.(varname)(:).*DA(:)));
    end
    
end

%% save results
info.td_vec = td_vec;
info.area = nansum(DA(:));
info.volume = nansum(DA(:).*G.h(:));

save([pth,'series_',island_tag,'.mat'],'P2','K2','E2','I2','SW','info');

