function [basename,nn_vec,dir0] = Z_runspec_raw
% 6/12/2014 Parker MacCready
%
% specify input-output locations, and file numbers to process
% 
% the intention is that you just change "calc_name" here when moving to
% different computers

calc_name = 'mac1';

switch calc_name
    case 'mac0' % single time, for testing
        basename = 'D2005';
        nn_vec = 1836; %
        dir1 = '/Users/PM5/Documents/roms/output/';
    case 'mac1' % a longer series, suitable for averaging
        basename = 'D2005';
        nn_vec = 1800:1928; % avg & dia - have 1800:1928 (his to 1929)
        dir1 = '/Users/PM5/Documents/roms/output/';
    case 'fjo' % fjord, for the full year
        basename = 'Cdia2005';
        nn_vec = 1:8760; % avg & dia - have 1:8760 (his to 8761)
        dir1 = '/data1/parker/roms/output/';
end

dir0.avg = [dir1,basename,'_avg/'];
dir0.dia = [dir1,basename,'_dia/'];
dir0.his = [dir1,basename,'_his/'];
