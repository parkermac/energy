function Z_res_plot(G,p2,k2,e2,sw,do_lp,Tdir,info)
% 6/13/2014  Parker MacCready

% figure
% Z_fig(14)
% set(gcf,'position',[250 10 2000 1300]);

subplot(121)
Z_pcolorcen(G.lon_rho,G.lat_rho,log10(p2.ape));
caxis([1 6]);
%colorbar('south')
title('log_{10} APE (J m^{-2})')
Z_dar;
Z_addcoast('combined',Tdir.coast);
[xt,yt] = Z_lab('ll');
text(xt,yt+1,datestr(info.td_avg))

subplot(122)
Z_pcolorcen(G.lon_rho,G.lat_rho,log10(k2.ke));
caxis([1 6]);
colorbar('south')
title('log_{10} KE (J m^{-2})')
Z_dar;
Z_addcoast('combined',Tdir.coast);
    
