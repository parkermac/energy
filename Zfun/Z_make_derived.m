function [p2,k2,e2,i2,sw,info] = Z_make_derived(p2,k2,sw,info)
% 5/21/2015  Parker MacCready
%
% derived quantities

% APE
p2.pdt = p2.rate;
p2.adv = p2.xadv + p2.yadv + p2.vadv;
p2.diss = p2.hdiff + p2.vdiff;
p2.err = p2.pdt - p2.adv - p2.diss - p2.dpe0dt;
% add the SW terms
p2.pdt_full = p2.pdt + sw.pdt;
p2.adv_full = p2.adv + sw.pdt;

% KE
k2.kdt = k2.accel;
k2.adv = k2.xadv + k2.yadv + k2.vadv + k2.cor + k2.prsgrd;
k2.diss = k2.vvisc;
k2.err = k2.kdt - k2.adv - k2.diss;

% KE + PE
e2 = struct();
e2.edt = p2.pdt_full + k2.kdt;
e2.bern = p2.adv_full + k2.adv;
e2.diss = p2.diss + k2.diss;
e2.background = p2.dpe0dt;
e2.err = e2.edt - e2.bern - e2.diss - e2.background;

% SW
sw.edt = sw.kdt + sw.pdt;
sw.diss = sw.sstr + sw.bstr;
sw.err = sw.edt - sw.bern - sw.diss;

% KE + PE - SW (baroclinic terms)
i2 = struct();
i2.kdt = k2.kdt - sw.kdt;
i2.pdt = p2.pdt_full - sw.pdt;
i2.edt = i2.kdt + i2.pdt;
i2.bern = e2.bern - sw.bern;
i2.diss = e2.diss - sw.diss;
i2.background = e2.background;
i2.err = e2.err - sw.err;
i2.wind_work = k2.wind_work;


