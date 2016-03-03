% Plot the integrating regions.

% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');
fn = '/Users/PM5/Documents/roms/output/D2005_his/ocean_his_1900.nc';

%% call the function
[G,S,T] = Z_get_basic_info(fn);
island0 = Z_island(G);

sal = load('Zfun/mask_salish.mat'); % has vectors x and y
insalish = inpolygon(G.lon_rho,G.lat_rho,sal.x,sal.y);


%% plotting

close all
figure
Z_fig(18)
set(gcf,'position',[10 10 1200 800]);
cax = [0 3000];
island_tag_list = {'full','shelf','salish'};
lab_list = {'(a)','(b)','(c)'};
for ii = 1:length(island_tag_list)
    island_tag = island_tag_list{ii};
    island = island0;
    % make a more restricted version of the mask
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
    hh = G.h;
    hh(island) = NaN;
    subplot(1,3,ii)
    cm = colormap('parula');
    cm = flipud(cm);
    colormap(cm)
    Z_pcolorcen(G.lon_rho,G.lat_rho,hh);
    caxis(cax);
    shading flat
    hold on; contour(G.lon_rho,G.lat_rho,G.h,[200 200],'-k')
    Z_dar;
    Z_addcoast('combined',Tdir.coast);
    title([lab_list{ii},' ',island_tag])
    xlabel('Longitude (deg)')
    if ii == 1
        ylabel('Latitide (deg)')
        [xt,yt] = Z_lab('lr');
        text(xt,yt,'Contour at 200 m','horizontalalignment','right')
    end
    
    if ii == 3
        colorbar('south');
        text(-126,44,'Bottom Depth (m)')
    end
        
end

%% Printing

if 1
    set(gcf,'PaperPositionMode','auto');
    print('-dpng',[Tdir.output,['energy_out/regions.png']]);
end




