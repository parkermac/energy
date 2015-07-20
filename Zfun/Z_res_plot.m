function Z_res_plot(G,p2,k2,e2,i2,sw,do_lp,Tdir,fn,info)
% 6/1/2015  Parker MacCready

colormap jet

p2.ape(p2.ape<0) = nan;
p2.ape_up(p2.ape_up<0) = nan;
p2.ape_down(p2.ape_down<0) = nan;
k2.ike =  k2.ke - sw.ke;
k2.ike(k2.ike<0) = nan;

NR = 1; NC = 3;

% make structures a and b, with fields and names, respectively
%
%a.ape = log10(p2.ape); b.ape = 'log_{10} iAPE (J m^{-2})';
a.ape_up = log10(p2.ape_up); b.ape_up = 'log_{10} iAPE-up (J m^{-2})';
a.ape_down = log10(p2.ape_down); b.ape_down = 'log_{10} iAPE-down (J m^{-2})';
%a.swape = log10(sw.pe); b.swape = 'log_{10} SW APE (J m^{-2})';
%a.ke = log10(k2.ke); b.ke = 'log_{10} KE (J m^{-2})';
a.ike = log10(k2.ike); b.ike = 'log_{10} iKE (J m^{-2})';
%a.swke = log10(sw.ke); b.swke = 'log_{10} SW KE (J m^{-2})';

fnm = fieldnames(a);
for ii = 1:length(fnm)
    varname = fnm{ii};
    subplot(NR,NC,ii)
    Z_pcolorcen(G.lon_rho,G.lat_rho,a.(varname));
    caxis([1 6]);
    hold on; contour(G.lon_rho,G.lat_rho,G.h,[200 200],'-k')
    if ii == length(fnm); colorbar('south'); end;
    title(b.(varname))
    Z_dar;
    Z_addcoast('combined',Tdir.coast);
    if ii == 1
        [xt,yt] = Z_lab('ll');
        text(xt,yt,datestr(info.td_avg))
    end
end
