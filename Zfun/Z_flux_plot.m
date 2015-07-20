function Z_flux_plot(G,p2,k2,e2,i2,sw,do_lp,Tdir,fn,info)
%function Z_res_plot(G,p2,k2,e2,i2,sw,do_lp,Tdir,info)
% 5/12/2015  Parker MacCready

colormap jet;
scl = .25;
if do_lp; scl = .05; end;
cax = scl*[-1 1];


sw_list = {'kdt','pdt','bern','diss'};
i2_list = {'kdt','pdt','bern','diss','wind_work'};

NR = 2;
NC = 6;

for ii = 1:length(sw_list)
    fld_name = sw_list{ii};
    fld = sw.(fld_name);
    subplot(NR,NC,ii)
    Z_pcolorcen(G.lon_rho,G.lat_rho,fld);
    caxis(cax);
    hold on; contour(G.lon_rho,G.lat_rho,G.h,[200 200],'-k')
    title(['sw ',strrep(fld_name,'_',' '),' (W m^{-2})'])
    Z_dar;
    Z_addcoast('combined',Tdir.coast);
end

for ii = 1:length(i2_list)
    fld_name = i2_list{ii};
    fld = i2.(fld_name);
    subplot(NR,NC,NC+ii)
    Z_pcolorcen(G.lon_rho,G.lat_rho,fld);
    caxis(cax);
    if ii == 1
        [xt,yt] = Z_lab('lr');
        if do_lp; tag = 'lp'; else; tag = 'raw'; end;
        ttext = strrep(fn,'.mat','');
        ttext = strrep(ttext,'flux_0','');
        text(xt,yt,[tag,' ',ttext], ...
            'horizontalalignment','r');
    end
    if ii == 4; colorbar('south'); end;
    hold on; contour(G.lon_rho,G.lat_rho,G.h,[200 200],'-k')
    title(['i2 ',strrep(fld_name,'_',' '),' (W m^{-2})'])
    Z_dar;
    Z_addcoast('combined',Tdir.coast);
end


