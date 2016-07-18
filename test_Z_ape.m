% Code to test the function Z_ape.

% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');
fn = '/Users/PM5/Documents/roms/output/D2005_his/ocean_his_4812.nc';

%% call the function
[G,S,T] = Z_get_basic_info(fn);

% test effect of region used for calculation
for nn = 10:10:110
    
    island = Z_island_nn(G, nn);
    [H] = Z_hyp(G,S,island);
    
    % Here is the function call.
    [apea, apea_up, apea_down] = Z_ape(fn,G,S,H);
    
    disp(['nn = ',num2str(nn),'  APEmax = ',num2str(round(nanmax(apea(:))/1000))])
    
end

%% plotting

if 0
    %close all
    figure
    Z_fig(14)
    set(gcf,'position',[10 10 1200 800]);
    cax = [2 6];
    fld_list = {'apea','apea_up','apea_down'};
    for ii = 1:length(fld_list)
        fld_name = fld_list{ii};
        eval(['fld = ',fld_name,';']);
        subplot(1,3,ii)
        colormap jet
        Z_pcolorcen(G.lon_rho,G.lat_rho,log10(fld));
        caxis(cax);
        shading flat
        colorbar
        hold on; contour(G.lon_rho,G.lat_rho,G.h,[200 200],'-k')
        title(['log10 ',strrep(fld_name,'_',' '),' (J m^{-2})'])
        Z_dar;
        Z_addcoast('combined',Tdir.coast);
    end
end


