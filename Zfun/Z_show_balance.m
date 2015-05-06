function Z_show_balance(p2,k2,e2,sw,info,DA)
% 6/12/2014  Parker MacCready

% do area integrals to get volume integrals (Watts)
% (or Joules for ape_avg)
P2 = struct();
K2 = struct();
E2 = struct();
SW = struct();

for vv = 1:length(info.ts_list)
    vname = info.ts_list{vv};
    P2.(vname) = squeeze(nansum(p2.(vname)(:).*DA(:)));
end
for vv = 1:length(info.uv_list)
    vname = info.uv_list{vv};
    K2.(vname) = squeeze(nansum(k2.(vname)(:).*DA(:)));
end
for vv = 1:length(info.e2_list)
    vname = info.e2_list{vv};
    E2.(vname) = squeeze(nansum(e2.(vname)(:).*DA(:)));
end
for vv = 1:length(info.sw_list)
    vname = info.sw_list{vv};
    SW.(vname) = squeeze(nansum(sw.(vname)(:).*DA(:)));
end

% screen output

scale = 1e6;
fmt = '%0.5g';
disp(['********************************'])
disp('P2')
disp([' dAPE/dt = ',sprintf(fmt,P2.rate/scale),' MW']);
disp([' check   = ',sprintf(fmt,P2.rate_check/scale),' MW']);
disp([' adv = ',sprintf(fmt,P2.adv/scale),' MW']);
disp([' diss = ',sprintf(fmt,P2.diss/scale),' MW']);
disp([' err = ',sprintf(fmt,P2.err/scale),' MW']);
disp('K2')
disp([' dKE/dt = ',sprintf(fmt,K2.accel/scale),' MW']);
disp(['check = ',sprintf(fmt,K2.accel_check/scale),' MW']);
disp([' adv = ',sprintf(fmt,K2.adv/scale),' MW']);
disp([' diss = ',sprintf(fmt,K2.diss/scale),' MW']);
disp([' err = ',sprintf(fmt,K2.err/scale),' MW']);
disp('SW')
disp(['d(APE+KE)/dt = ',sprintf(fmt,SW.aedt/scale),' MW']);
disp([' check       = ',sprintf(fmt,SW.aedt_check/scale),' MW']);
disp([' abern = ',sprintf(fmt,SW.abern/scale),' MW']);
disp([' diss = ',sprintf(fmt,SW.diss/scale),' MW']);
disp([' aerr = ',sprintf(fmt,SW.aerr/scale),' MW']);


