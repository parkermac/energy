% set up the environment
clear; addpath('../alpha/'); Tdir = toolstart;
addpath('./Zfun');


[basename,nn_vec,dir0] = Z_runspec_raw; % just to get dir0
nn = 4800;
ns = num2str(nn); ns = ['0000',ns]; ns = ns(end-3:end);
fn_avg = [dir0.avg,'ocean_avg_',ns,'.nc'];

odir = [Tdir.output,'energy_out/'];

load([odir,'mix_out_',num2str(nn),'_71.mat']);

%% Get hypsometry (only needed for the island mask)
[G,S,T] = Z_get_basic_info(fn_avg);
island = Z_island(G);
[H] = Z_hyp(G,S,island);
[zr,zw] = Z_s2z(G.h,0*G.h,S);

mix3(H.island3) = NaN;

zr(H.island3) = NaN;
G.h(H.island) = NaN;
for kk = 1:S.N+1
    this_zw = zw(kk,:,:);
    this_zw(H.island) = NaN;
    zw(kk,:,:) = this_zw;
end


%% plotting
close all
figure
Z_fig(16)
set(gcf,'position',[10 10 1300 1000]);

%% map

lat_list = [49.2, 47, 44.3];

ax_mat = [-124 -123 -400 0.1; -125 -124 -100 0.1; -125 -124 -100 0.1];

subplot(121);
colormap jet
Z_pcolorcen(G.lon_rho,G.lat_rho,1e3*mix2);
caxis(10*[-1 1]);
shading flat
colorbar('south')
hold on; contour(G.lon_rho,G.lat_rho,G.h,[200 200],'-k')
title('(a) Mixing_{A} (mW m^{-2})')
Z_dar;
Z_addcoast('combined',Tdir.coast);

xlabel('Longitude (deg)')
ylabel('Latitude (deg)')

lab_list = ['b','c','d'];
for ii = 1:length(lat_list)
    lat = lat_list(ii);
    jlat = find(G.lat_rho(:,1) >= lat); jlat = jlat(1);
    llat = G.lat_rho(jlat,1);
    hold on
    plot([ax_mat(ii,1) ax_mat(ii,2)],[llat llat], '-k','linewidth',3)
    text(ax_mat(ii,1)-.2,llat,lab_list(ii))
end

%% sections

if 1
    
    for ii = 1:length(lat_list)
        lat = lat_list(ii);
        jlat = find(G.lat_rho(:,1) >= lat); jlat = jlat(1);
        Lat = G.lat_rho(jlat,:);
        Lon = G.lon_rho(jlat,:);
        Mix3 = squeeze(mix3(:,jlat,:));
        Sig = squeeze(den1(:,jlat,:));
        Zr = squeeze(zr(:,jlat,:));
        Zw = squeeze(zw(:,jlat,:));
        
        Zf = [Zw(1,:) ; Zr ; Zw(end,:)];
        Mix3f = [Mix3(1,:) ; Mix3 ; Mix3(end,:)];
        Sigf = [Sig(1,:) ; Sig ; Sig(end,:)];
        Lonf = ones(size(Zf,1),1) * Lon;
        
        subplot(length(lat_list),2,2*ii)
        pcolor(Lonf,Zf,1e3*Mix3f)
        caxis(1*[-1 1])
        shading interp
        colorbar('eastoutside')
        
        hold on
        contour(Lonf, Zf, Sigf,[20:.5:30],'-k')
        
        plot(Lon,Zw(1,:),'-k','linewidth',3)
        
        axis(ax_mat(ii,:))
        
        xlabel('Longitude (deg)')
        ylabel('Z (m)')
        if ii == 1
            title('(b) Mixing_{V} (mW m^{-3})')
        elseif ii == 2
            [xt,yt] = Z_lab('ll');
            text(xt,yt,'(c)','fontweight','bold','color','w')
        elseif ii == 3
            [xt,yt] = Z_lab('ll');
            text(xt,yt,'(d)','fontweight','bold','color','w')
            
        end
        
        box on
    end
    
end

%% printing

if 0
    set(gcf,'PaperPositionMode','auto');
    print('-dpng',[Tdir.output,['energy_out/mix_sections.png']]);
end


