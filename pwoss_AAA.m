see_filters = 0;
load_data = 1;
if load_data == 1  
    load rfc_big_scan_170927.mat;
    load rfc_FL_scan_170927.mat;
    load cones55555.mat;
    load kLOS_171122.mat;
    load bcb_findss.mat;
end
pix2eV = (1.599/(2*496));
pix2invA = 0.512*0.04631/180*3.1415*14/30*sqrt(110-4);

AVE_FL = mean(rfc_FL_Es(rfc_FL_Es>0));


table_title = 'BCB fit';
dirac_Es = bcb_finds_180109;%502*ones(1,961);% BCB_Es;%BCB_Es;%rfc_Es_after;   %input Energies vector (cone pixels)
dirac_ks = dos_Ks_180109;%181*ones(1,961);%kLOS;%181*ones(1,961);% kLOS;%.5*(rfc_ks_after+kLOS);   %input mtm vector (cone pixels)

DPEs_eV = (AVE_FL - dirac_Es)*pix2eV;%(AVE_FL - dos_Es)*pix2eV;%(rfc_FL_Es - bcb_findss)*pix2eV;%(645*ones(1,961)-dirac_Es)*pix2eV;%(rfc_FL_Es - bcb_findss)*pix2eV;%bcb_finds;%rfc_Es_after;%(rfc_FL_Es - rfc_small_cut_Es)*pix2eV;
B_map = reshape(DPEs_eV,31,31);   

B_interval_list_abs =[.14,.16; .16,.18; .18,.20; .20,.22; .22,.24; .24,.26; .26,.28];%[.15,.18; .18,.2; .2,.22; .22,.25; .25,.28];% [.462,.4674; .4674,.4701; .4701,.4728; .4728,.4782; .4782,.4836; .4836,.489];%[.2,.21; .21,.22; .22, .23; .24,.25; .25,.26; .26,.27];%[.22,.225; .225,.23; .23,.235; .235,.24; .24,.25; .25,.26];%[.22,.2236; .2236,.226; .2271,.2283; .2283,.233];%[400,500; 500,600];%[0,.5; .5,1];%.195,.205; .205,.215; .215,.225; .225,.235; .235,245];%; 480,500; 500,510; 510,520; 520,535; 535,555];% 495,500; 505,510; 510,525; 525,540; 540,555];%[320,340; 340,350; 350,355; 355,360; 360,380];
[.15,.18; .18,.2; .2,.22; .22,.25; .25,.28];
B_interval_list = (B_interval_list_abs - min(B_map(:)))/(max(B_map(:))-min(B_map(:)));

figure()
E_hist = histogram(DPEs_eV(DPEs_eV<0.9),40); hold on;
B_hist_list = reshape(B_interval_list_abs,1,[]);
plot([B_hist_list],[ones(size(B_hist_list))],'r*'), hold off;


pre_filter = cat(1,[]);
pre_filter(1,:) = [172,188,kLOS];%[170,189, kLOS];
pre_filter(2,:) = [220,440, DPI_big];

allfiltered = ones(1,961);
for NN = 1:size(pre_filter,1)
    allfiltered(pre_filter(NN,3:end) < pre_filter(NN,1)) = 0;
    allfiltered(pre_filter(NN,3:end) > pre_filter(NN,2)) = 0;
    if see_filters == 1
        if NN == 1
            figure
        end
        subplot(size(pre_filter,1),2,2*(NN-1)+2)
        histogram(pre_filter(NN,3:end),40), hold on;
        plot([pre_filter(NN,1),pre_filter(NN,1)],[0,max(hist(pre_filter(NN,3:end),40))],'r'), hold on;
        plot([pre_filter(NN,2),pre_filter(NN,2)],[0,max(hist(pre_filter(NN,3:end),40))],'r'), hold off;            
    end
    if NN == size(pre_filter,1)
        filtered_A_map = reshape(allfiltered,31,31);
        if see_filters==1
            subplot(NN,2,[1]) 
            imagesc(filtered_A_map), axis xy
        end
        A_map = filtered_A_map;
        A_interval_list = [0.5,1];
    end
end


[E_interp,K_interp] = meshgrid(-399:1:350, -100:1:100); 
full_arpes = zeros(size(E_interp));

[region_list]=Correlation_list_indep_ranges(A_map,B_map,A_interval_list,B_interval_list);
regional_arpes=zeros([size(E_interp),size(region_list,1)]);

involved_scans = zeros(1,961);
involved_scans_panel = zeros(1,961,size(region_list,1));
mean_DPE_eV = zeros(1,size(region_list,1));
mean_FL_pix = zeros(1,size(region_list,1));
mean_FL_eV = zeros(1,size(region_list,1));
map_regions = cat(3,[]);
NNN = 0;
for iii=1:size(region_list,1)
    nnn = 0;
    
     for jjj= find(region_list{iii}==1)'
         involved_scans(jjj) = involved_scans(jjj)+1;
         involved_scans_panel(1,jjj,iii) = involved_scans_panel(1,jjj,iii)+1;
         nnn = nnn+1;
         [map_y, map_x] = ind2sub([31,31],jjj);
         map_regions(1,nnn,iii) = map_x;
         map_regions(2,nnn,iii) = map_y;
         
         cone = result1i(:,:,jjj);
         [E_coor, K_coor] = meshgrid(1:size(cone,2),1:size(cone,1));
         E_coor = E_coor - dirac_Es(jjj);  
         K_coor = K_coor - dirac_ks(jjj);  
         interp_arpes = interp2(E_coor, K_coor, cone, E_interp, K_interp);
         figure, imagesc(interp_arpes), axis xy
         regional_arpes(:,:,iii) = regional_arpes(:,:,iii) + interp_arpes;
         
         DPE_pix = rfc_FL_Es(jjj);
         FL_pix = rfc_FL_Es(jjj);
         DPE_eV = (-FL_pix + DPE_pix) * pix2eV;
         
         DPE_eV = (AVE_FL - dirac_Es(jjj))*pix2eV;
         FL_eV = FL_pix * pix2eV;
         mean_DPE_eV(iii) = mean_DPE_eV(iii) + DPE_eV;
         mean_FL_pix(iii) = mean_FL_pix(iii) + FL_pix;
         mean_FL_eV(iii) = mean_FL_eV(iii) + FL_eV;
     end
  
          
     mean_DPE_eV(iii) = mean_DPE_eV(iii) / length(find(region_list{iii}==1));
     mean_FL_pix(iii) = mean_FL_pix(iii) / length(find(region_list{iii}==1));
     mean_FL_eV(iii) = mean_FL_eV(iii) / length(find(region_list{iii}==1));
     nnn_scans(iii) = nnn;
     NNN = NNN + nnn;
end


plot out the binned and symmetrized regional spectra

E_bin=5;
K_bin=3;

E_binned = E_interp(1,1:E_bin:end-1);

if K_bin == 1
    K_binned = K_interp(:,1);
else    
    K_binned = K_interp(1:K_bin:end,1);
end

Eaxis=E_binned;
arp_chart = figure;
nn=0;
nnn=0;
out_spec = cell([1,size(region_list,1)]);
out_kax = cell([1,size(region_list,1)]);
out_eax = cell([1,size(region_list,1)]);
out_nnn = zeros(1,size(region_list,1));
for iii=1:size(region_list,1)    
    regional_arpes(:,:,iii) = imgaussfilt(regional_arpes(:,:,iii),2.5);    
    region_arpes_binned = Binning_2d(regional_arpes(:,:,iii),E_bin,K_bin);
    [region_arpes_symmetrized,Kaxis]=Symmetrized_spectra(region_arpes_binned,K_binned);
    ax = subplot(2,size(B_interval_list,1),iii); %%%%%%replace 2 with size(B_interval_list,1)%%%
    
    imagesc(regional_arpes(:,:,iii)), axis xy
    imagesc(rot90(Eaxis,-1),rot90(Kaxis,-1),rot90(region_arpes_binned(1:103,:),-1)), %axis xy;
    imagesc(Eaxis,Kaxis,region_arpes_symmetrized), axis xy;
    Eaxis_eV = E_interp(1,1:E_bin:end-1) * pix2eV;
    Eaxis_eV_0fixed = Eaxis_eV - abs(mean_DPE_eV(iii));           
    Kaxis_invA = Kaxis * pix2invA;
    size(Kaxis_invA)
    
    if length(find(region_list{iii}==1)) > 0
        imagesc(Kaxis_invA, flip(Eaxis_eV_0fixed), rot90(norman(region_arpes_symmetrized,0,5))), axis xy
        imagesc(Kaxis_invA, flip(Eaxis_eV_0fixed), rot90(norman(region_arpes_binned(1:103,:),0,5))), axis xy
        imagesc(Kaxis_invA, flip(Eaxis_eV_0fixed), rot90(norman(region_arpes_binned,0,5))), axis xy
        imagesc(Kaxis_invA,flip(Eaxis_eV_0fixed), rot90(norman(regional_arpes(:,:,iii),0,5))), axis xy;
        
        ylim([0,700])
        ylim([-.8,0.25])
        yticks([0,-(mean_DPE_eV(iii)/pix2eV)])
        yticklabels({num2str(mean_DPE_eV(iii)),'0'})
        
        yt = get(gca, 'YTick');
        set(gca, 'FontSize',7)
        ytickangle(90)
        
        xticks([-50,0,50])
        xticklabels({num2str(-50*pix2invA),'0',num2str(50*pix2invA)})
        xlabel('Momentum ({\AA}$^{-1}$)','Interpreter','latex','FontSize',8)
        xt = get(gca, 'XTick');
        set(gca, 'FontSize', 7)
    end
    
    if iii == 1+size(B_interval_list,1)*nn
        descrr = {[num2str(A_interval_list(nn+1,1)),'-',num2str(A_interval_list(nn+1,2))]};
        axes(ax)
        ylabel_str = 'E - E$_{F}$, (eV)';
        ylabel(ylabel_str,'Interpreter','latex','FontSize',8)
        text(-0.2,-0.3,descrr')
        ylabel([num2str(A_interval_list(nn+1,1)),'-',num2str(A_interval_list(nn+1,2))]);
        nn=nn+1;
    end 
    
    if ismember(iii,[1:size(B_interval_list,1)]) == 1
        set(gca,'XaxisLocation','top');
        descr = {[num2str(B_interval_list_abs(nnn+1,1)),'-',num2str(B_interval_list_abs(nnn+1,2))]};
        axes(ax)
        text(-.05,800,descr)
        text(-.05,0.35,descr)
        xlabel([num2str(B_interval_list(nnn+1,1)),'-',num2str(B_interval_list(nnn+1,2))]);
        nnn=nnn+1;
    end
    
    ax2 = subplot(2,size(B_interval_list,1),iii+size(B_interval_list,1));
    imagesc(B_map), axis xy, hold on; %caxis([min(B_map(B_map>0)),max(B_map(B_map>0))]), hold on;
    plot(map_regions(1,:,iii),map_regions(2,:,iii),'w+'), hold off;
    title(['nnn=',num2str(nnn_scans(iii))])
    colormap jet
    
    output_kaxis = [1:size(region_arpes_binned,1)] - round(size(region_arpes_binned,1)/2);
    output_kaxis_invA = output_kaxis * K_bin * pix2invA;
    
    out_spec{1,iii} = regional_arpes(:,:,iii);
    out_kax{1,iii} = output_kaxis_invA;
    out_eax{1,iii} = Eaxis_eV_0fixed;
    out_nnn(iii) = nnn_scans(iii);
end
suptitle([table_title,', total scans: ',num2str(NNN)]);





