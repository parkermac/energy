% sw_fix.m  6/18/2014  Parker MacCready
%
% patches over (7) jumps in the sw flux fields

clear
addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');

pth = '/Users/PM5/Documents/tools_output/energy_out/Cdia2005/flux_raw/';
nn_vec = 1:8760;
NT = length(nn_vec);

load([pth,'G.mat']); % has island

tt = 0; % a counter
for nn = nn_vec
    
    tt = tt + 1; % increment counter
    
    ncpad = ['0000',num2str(nn)]; ncpad = ncpad(end-4:end);
    infile = ['flux_',ncpad,'.mat'];
    %if(mod(nn,10)==0); disp(['Working on ',infile]); end;
    load([pth,infile]);
    
    if tt==1
        G.island = isnan(k2.accel_check); % redefine G.island
        DA = G.DX .* G.DY;
        DA(G.island) = NaN;
    end
    
    sw_aedt = sw.kdt + sw.apdt;
    sw_diss = sw.sstr + sw.bstr;
    sw_aerr = sw_aedt - sw.abern - sw_diss;
    SW_aerr = squeeze(nansum(sw_aerr(:).*DA(:)));
    
    if abs(SW_aerr) > 1e4;
        
        disp(['Patching file ',num2str(nn)]);
        
        % save a copy of the old file
        movefile([pth,infile],[pth,'ORIG_',infile]);
        
        % name the before and after files
        ncpad1 = ['0000',num2str(nn-1)]; ncpad1 = ncpad1(end-4:end);
        infile1 = ['flux_',ncpad1,'.mat'];
        ncpad2 = ['0000',num2str(nn+1)]; ncpad2 = ncpad2(end-4:end);
        infile2 = ['flux_',ncpad2,'.mat'];
        
        if tt==1
            f2 = load([pth,infile2]);
        elseif tt==NT;
            f1 = load([pth,infile1]);
        else
            f1 = load([pth,infile1]);
            f2 = load([pth,infile2]);
        end

        for vv = 1:length(info.sw_list)
            vname = info.sw_list{vv};
            if tt==1
                sw.(vname) = f2.sw.(vname);
            elseif tt==NT;
                sw.(vname) = f1.sw.(vname);
            else
                sw.(vname) = (f1.sw.(vname) + f2.sw.(vname))/2;
            end
        end

        save([pth,infile],'p2','k2','sw','info');
    end
end


