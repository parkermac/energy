%% flux_res_plot.m  5/12/2015  Parker MacCready
%
% plots fluxes or reservoirs

clear
addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');

do_flux = 1;

if 1 % user selection
    pth0 = '/Users/PM5/Documents/tools_output/energy_out/Cdia2005/';
    [fn,pth] = uigetfile(pth0,'Select .mat file(s) to plot', ...
        'multiselect','on');
else % for plot development
    fn = 'flux_01836.mat';
    pth = ['/Users/PM5/Documents/tools_output/energy_out/', ...
        'Cdia2005/flux_lp71/']; % _raw or _lp71
end

lp_test = strfind(pth,'lp71');
if ~isempty(lp_test)
    do_lp = 1;
else
    do_lp = 0;
end

load([pth,'G.mat']); % has island

make_movie = 0;
if iscell(fn); make_movie = 1; end
% determine how many files to plot
if make_movie; ntt = size(fn,2); else ntt = 1; end;

figure
Z_fig(14);
set(gcf,'position',[1 1 1400 825]); % laptop version
%set(gcf,'position',[10 10 2500 1300]);
set(gcf,'PaperPositionMode','auto');

%%
for tt = 1:ntt % MOVIE loop start (or just make single plot)
    if make_movie;
        fnt = fn{tt};
        infile = [pth,fnt];
        outfile = strrep(fnt,'.mat','.jpg');
    else
        fnt = fn;
        infile = [pth,fnt];
        outfile = strrep(fnt,'.mat','.jpg');
    end
    load(infile);
    
    [p2,k2,e2,i2,sw,info] = Z_make_derived(p2,k2,sw,info);
    
    if do_flux
        Z_flux_plot(G,p2,k2,e2,i2,sw,do_lp,Tdir,fnt);
        odir_name = 'a_flux_movie';
    else
        Z_res_plot(G,p2,k2,e2,sw,do_lp,Tdir,info);
        odir_name = 'a_res_movie';
    end
    
    if make_movie == 0
        if 1
            % do nothing
        else
            % choose to save a copy of the plot
            plotit = input('Save a plot? [1 = save, RETURN = do not save]');
            if isempty(plotit); plotit = 0; end;
            if plotit
                print('-djpeg100',[fn,'.jpg']);
            end
        end
    else % make a folder of jpegs for a movie
        if tt==1
            outdir = [pth,odir_name,'/'];
            if exist(outdir,'dir')==7; rmdir(outdir,'s'); end;
            mkdir(outdir);
        end
        print('-djpeg100',[outdir,outfile]);
        if tt<length(fn); clf; end;
    end
    
end
 