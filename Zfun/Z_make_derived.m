function [p2,k2,e2,sw,info] = Z_make_derived(p2,k2,sw,info)
% 5/5/2015  Parker MacCready
%
% derived quantities

% PE
p2.adv = p2.xadv + p2.yadv + p2.vadv;
p2.diss = p2.hdiff + p2.vdiff;
p2.err = p2.rate - p2.adv - p2.diss - p2.dpe0dt;
info.ts_list{length(info.ts_list)+1} = 'adv';
info.ts_list{length(info.ts_list)+1} = 'diss';
info.ts_list{length(info.ts_list)+1} = 'err';

% KE
k2.adv = k2.xadv + k2.yadv + k2.vadv + k2.cor + k2.prsgrd;
k2.diss = k2.vvisc;
k2.err = k2.accel - k2.adv - k2.diss;
info.uv_list{length(info.uv_list)+1} = 'adv';
info.uv_list{length(info.uv_list)+1} = 'diss';
info.uv_list{length(info.uv_list)+1} = 'err';

% KE + PE
e2 = struct();
e2.edt = p2.rate + k2.accel;
e2.bern = p2.adv + k2.adv;
e2.diss = p2.diss + k2.diss;
e2.background = p2.dpe0dt;
e2.err = e2.edt - e2.bern - e2.diss - e2.background;
info.e2_list = {'edt','bern','diss','background','err'};

% SW
sw.edt = sw.kdt + sw.pdt;
sw.edt_check = sw.kdt_check + sw.pdt_check;
sw.aedt = sw.kdt + sw.apdt;    
sw.aedt_check = sw.kdt_check + sw.apdt_check;
sw.diss = sw.sstr + sw.bstr;
sw.err = sw.edt - sw.bern - sw.diss;
sw.aerr = sw.aedt - sw.abern - sw.diss;
info.sw_list{length(info.sw_list)+1} = 'edt';
info.sw_list{length(info.sw_list)+1} = 'edt_check';
info.sw_list{length(info.sw_list)+1} = 'aedt';
info.sw_list{length(info.sw_list)+1} = 'aedt_check';
info.sw_list{length(info.sw_list)+1} = 'diss';
info.sw_list{length(info.sw_list)+1} = 'err';
info.sw_list{length(info.sw_list)+1} = 'aerr';
