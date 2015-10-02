% This is designed to test the calculation of dAPEa/dt
%
% file timing and numbering:
% nn     o   nn+1
% his    o   his
%  |    nn    |
%  |  avg,dia |

%% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');

dir0.his = '/Users/PM5/Documents/roms/output/D2005_his/';
dir0.avg = '/Users/PM5/Documents/roms/output/D2005_avg/';
dir0.dia = '/Users/PM5/Documents/roms/output/D2005_dia/';

for nn = 1800:1824
    ns = num2str(nn); ns = ['0000',ns]; ns = ns(end-3:end);
    fn = [dir0.his,'ocean_his_',ns,'.nc'];
    [G,S,T] = Z_get_basic_info(fn);
    island = Z_island(G);
    [H] = Z_hyp(G,S,island);
    % call functions
    p2 = new_flux(dir0,nn,H);
    p2.rate_error = p2.rate - p2.rate_check;
    disp('-----------------------------------')
    disp(['nn = ',num2str(nn)])
    disp(['std(dAPE/dt) = ',num2str(nanstd(p2.rate(:))),' W'])
    disp(['std(Error)   = ',num2str(nanstd(p2.rate_error(:))),' W'])
end
%% plotting

if 0   
    %close all
    figure
    Z_fig(14)
    set(gcf,'position',[10 10 1200 800]);
    
    colormap jet;
    scl = .5;
    cax = scl*[-1 1];
    
    fld_list = {'rate', 'rate_check', 'rate_error','dapedt_uno','dapedt_dos'};
    for ii = 1:length(fld_list)
        fld_name = fld_list{ii};
        eval(['fld = p2.',fld_name,';']);
        subplot(2,3,ii)
        Z_pcolorcen(G.lon_rho,G.lat_rho,fld);
        caxis(cax);
        shading flat
        colorbar
        hold on; contour(G.lon_rho,G.lat_rho,G.h,[200 200],'-k')
        title([strrep(fld_name,'_',' '),' (W m^{-2})'])
        Z_dar;
        Z_addcoast('combined',Tdir.coast);
    end
end


