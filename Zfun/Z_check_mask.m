function Z_check_mask(G,info,p2,k2,sw,ISN,island_orig,Tdir)
% 6/13/2014 Parker MacCready

disp('ISLAND MASK')
disp([num2str(ISN),' masked points']);
disp('*** PE ***')
for vv = 1:length(info.ts_list)
    vname = info.ts_list{vv};
    isn = sum(isnan(p2.(vname)(:)));
    disp([num2str(isn),' NaNs in p2.',vname, ...
        ' (',num2str(isn-ISN),' extras)']);
end
disp('*** KE ***')
for vv = 1:length(info.uv_list)
    vname = info.uv_list{vv};
    isn = sum(isnan(k2.(vname)(:)));
    disp([num2str(isn),' NaNs in k2.',vname, ...
        ' (',num2str(isn-ISN),' extras)']);
end
disp('*** SW ***')
for vv = 1:length(info.sw_list)
    vname = info.sw_list{vv};
    isn = sum(isnan(sw.(vname)(:)));
    disp([num2str(isn),' NaNs in sw.',vname, ...
        ' (',num2str(isn-ISN),' extras)']);
end

% plotting

figure
Z_fig(14);

bvn = 'accel_check';
bv = k2.(bvn);
bad = isnan(bv) & ~island_orig;
[rbad,cbad] = find(bad);

Z_pcolorcen(G.lon_rho,G.lat_rho,bv);
caxis([-10 10]);
hold on
plot(G.lon_rho(bad),G.lat_rho(bad),'*r')
plot(G.lon_rho(bad),G.lat_rho(bad),'ok')
title([strrep(bvn,'_',' '),' mask errors'])
Z_dar;
Z_addcoast('combined',Tdir.coast);


