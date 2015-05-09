function Z_flux_plot(G,p2,k2,e2,sw,do_lp,Tdir)
% 5/5/2015  Parker MacCready

colormap jet;
scl = 1;
if do_lp; scl = .05; end;
cax = scl*[-1 1];

e2_list = {'edt','bern','diss','background'};

i2.edt = e2.edt - sw.aedt;
i2.bern = e2.bern - sw.abern;
i2.diss = e2.diss - sw.diss;
i2.background = e2.background;
i2.err = e2.err - sw.err;

for ii = 1:length(e2_list)
    fld_name = e2_list{ii};
    fld = i2.(fld_name);
    subplot(1,4,ii)
    Z_pcolorcen(G.lon_rho,G.lat_rho,fld);
    caxis(cax);
    if ii == 1; colorbar('south'); end;
    title([strrep(fld_name,'_',' '),' (W m^{-2})'])
    Z_dar;
    Z_addcoast('combined',Tdir.coast);

end

