%% flux_lp.m  6/12/2014  Parker MacCready
%
% does low-pass filtering of the vertically-integrated flux terms
% calculated by flux.m (via flux_driver.m)

clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');
odir_top = [Tdir.output,'energy_out/'];

% set indices of times to average around
nn_center_vec = 1848;

% set size of averaging window
win_len = 71;

% input directory
[basename,nn_vec,dir0] = Z_runspec_raw;
idir = [odir_top,basename,'/flux_raw/'];

% create EMPTY result directory, specific to this window length
odir = [odir_top,basename,'/flux_lp',num2str(win_len),'/'];
if exist(odir,'dir'); rmdir(odir,'s'); end; mkdir(odir);
disp(' '); disp(['Results in: ',odir]);

% get the grid (including the island mask)
load([idir,'G.mat']);
DA = G.DX .* G.DY;
DA(G.island) = NaN;

save([odir,'G.mat'],'G');

%%

% create weighting vector
if win_len == 71; wt = Z_godin_shape;
else; wt = hanning(win_len); wt = wt/sum(wt); end;

for cc = nn_center_vec % start of "day" loop
    
    tic
    
    ccpad = ['0000',num2str(cc)]; ccpad = ccpad(end-4:end);
    fn_out = [odir,'flux_',ccpad,'.mat'];
    
    nn1 = cc - floor(win_len/2);
    nn2 = nn1 + win_len - 1;
    disp([' nn1=',num2str(nn1),' nn2=',num2str(nn2)]);
    nn_vec = nn1:nn2;
    
    ww = 0; % a counter
    for nn = nn_vec % start of "hour" loop
        
        ww = ww+1;
        
        % name the infile
        nnpad = ['0000',num2str(nn)]; nnpad = nnpad(end-4:end);
        fn_in = [idir,'flux_',nnpad,'.mat'];
        
        if ww==1
            P2 = struct();
            K2 = struct();
            SW = struct();
        end
                    
        load(fn_in)
        
        % loop over the variables
        for ii = 1:length(info.ts_list)
            var = info.ts_list{ii};
            if ww==1
                P2.(var) = wt(ww)*p2.(var);
            else
                P2.(var) = P2.(var) + wt(ww)*p2.(var);
            end
        end
        for ii = 1:length(info.uv_list)
            var = info.uv_list{ii};
            if ww==1
                K2.(var) = wt(ww)*k2.(var);
            else
                K2.(var) = K2.(var) + wt(ww)*k2.(var);
            end
        end
        for ii = 1:length(info.sw_list)
            var = info.sw_list{ii};
            if ww==1
                SW.(var) = wt(ww)*sw.(var);
            else
                SW.(var) = SW.(var) + wt(ww)*sw.(var);
            end
        end
     
    end
    
    p2 = P2;
    k2 = K2;
    sw = SW;
    
    % note: we should fix this so it saves info for the center file.
    % so that the time stamp matches the center, instead of the end
    save(fn_out,'p2','k2','sw','info')
        
    dt = toc; disp(['   - took ',num2str(round(dt)),' sec'])
    
end % end of day loop (dd)








