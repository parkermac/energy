% Make the tidally-averaged the Mixing field
% and save 3D to show depth structure.

% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');

% set index of times to average around
nn_center_vec = 4620;

% set size of averaging window (71)
win_len = 71;

[basename,nn_vec,dir0] = Z_runspec_raw; % just to get dir0
nn = nn_center_vec(1);
ns = num2str(nn); ns = ['0000',ns]; ns = ns(end-3:end);
fn_avg = [dir0.avg,'ocean_avg_',ns,'.nc'];

odir = [Tdir.output,'energy_out/'];


%% Get hypsometry (only needed for the island mask)
[G,S,T] = Z_get_basic_info(fn_avg);
island = Z_island(G);
[H] = Z_hyp(G,S,island);
[zr,zw] = Z_s2z(G.h,0*G.h,S);

%% time averaging

% create weighting vector
if win_len == 71; wt = Z_godin_shape;
else; wt = hanning(win_len); wt = wt/sum(wt); end;

for cc = nn_center_vec; % start of "day" loop
    
    tic
    
    nn1 = cc - floor(win_len/2);
    nn2 = nn1 + win_len - 1;
    disp([' nn1=',num2str(nn1),' nn2=',num2str(nn2)]);
    NN_vec = nn1:nn2;
    
    ww = 0; % a counter
    for nn = NN_vec % start of "hour" loop
        
        ww = ww+1;
        
        if ww==1
            P2 = struct();
            P3 = struct();
        end
        
        % get the 3D APE fields (just mixing)
        [p2,p3,den1,info] = flux3d(dir0,nn,H,zr,zw);
        
        if nn == cc % save the time of the center of the window
            center_info = info;
        end
        
        % loop over the variables
        fnm = fieldnames(p2);
        for ii = 1:length(fnm)
            varname = fnm{ii};
            if ww==1
                P2.(varname) = wt(ww)*p2.(varname);
            else
                P2.(varname) = P2.(varname) + wt(ww)*p2.(varname);
            end
        end
        % loop over the variables
        fnm = fieldnames(p3);
        for ii = 1:length(fnm)
            varname = fnm{ii};
            if ww==1
                P3.(varname) = wt(ww)*p3.(varname);
            else
                P3.(varname) = P3.(varname) + wt(ww)*p3.(varname);
            end
        end
        % single field - no loop needed
        if ww==1
            DEN1 = wt(ww)*den1;
        else
            DEN1 = DEN1 + wt(ww)*den1;
        end
     
    end
    
    p2 = P2;
    p3 = P3;
    den1 = DEN1;
    
    info = center_info;
            
    dt = toc; disp(['   - took ',num2str(round(dt)),' sec'])
    
end % end of day loop (dd)


%% make mixing
mix2 = p2.vdiff + p2.hdiff;
mix3 = p3.vdiff + p3.hdiff;

%% save output for later plotting

save([odir,'mix_out_',num2str(nn_center_vec),'_',num2str(win_len),'.mat'],'mix2', 'mix3', 'den1');


