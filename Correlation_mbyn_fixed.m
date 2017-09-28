table_title = 'rfc small scan 170928';
dirac_Es = rfc_small_Es;   %input Energies vector (cone pixels)
dirac_ks = rfc_small_ks;   %input mtm vector (cone pixels)

B_map = E_small_map;      
B_interval_list = [.25,.3; .3,.35; .4,.45; .45,.5; .5,.55; .55, .6];% .275,.3; .3,.325; .325,.35]+.1;

pre_filter = cat(1,[]);
pre_filter(1,:) = [220,276, scan_spread];
pre_filter(2,:) = [172,188, rfc_small_ks];
pre_filter(3,:) = [168,434, DPI_big];
pre_filter(4,:) = [.85,1, scan_maxes_FL];

figure
allfiltered = reshape(B_map,1,961);
for NN = 1:size(pre_filter,1)
    allfiltered(pre_filter(NN,3:end) < pre_filter(NN,1)) = 0;
    allfiltered(pre_filter(NN,3:end) > pre_filter(NN,2)) = 0;
    
    subplot(size(pre_filter,1),2,2*(NN-1)+2)
    histogram(pre_filter(NN,3:end),40), hold on;
    plot([pre_filter(NN,1),pre_filter(NN,1)],[0,max(hist(pre_filter(NN,3:end),40))],'r'), hold on;
    plot([pre_filter(NN,2),pre_filter(NN,2)],[0,max(hist(pre_filter(NN,3:end),40))],'r'), hold off;
    
end

filtered_A_map = reshape(allfiltered,31,31);
subplot(NN,2,[1 3 5]) 
imagesc(filtered_A_map), axis xy
 

A_interval_list = [.5,1];
A_map = filtered_A_map;%eshape(norman(draw_rats_908,0,5),31,31);%reshape(line_ks_829,31,31); %reshape(ABEK_ks,31,31);%DPI_map;   %Horizontal axis




pix2eV = (1.599/(2*496));
pix2invA = 0.512*0.04631/180*3.1415*14/30*sqrt(110-4);


[E_interp,K_interp] = meshgrid(-149.5:.5:350, -80:.5:80); %meshgrid(-149.5:.5:280, -80:.5:80); %meshgrid(-199.5:0.5:250,-80:0.5:80);
full_arpes = zeros(size(E_interp));

[region_list]=Correlation_list_indep_ranges(A_map,B_map,A_interval_list,B_interval_list);
regional_arpes=zeros([size(E_interp),size(region_list,1)]);
%E_coor_ave = zeros(size(cone,2),size(region_list,1));
%K_coor_ave = zeros(size(cone,1),size(region_list,1));

mean_DPE_eV = zeros(1,size(region_list,1));
mean_FL_eV = zeros(1,size(region_list,1));
map_regions = cat(3,[]);
NNN = 0;
for iii=1:size(region_list,1)
    nnn = 0;
    
     for jjj=find(region_list{iii}==1)'
         
         nnn = nnn+1;
         [map_y, map_x] = ind2sub([31,31],jjj);
         map_regions(1,nnn,iii) = map_x;
         map_regions(2,nnn,iii) = map_y;
         cone = cones(:,:,jjj);
         [E_coor, K_coor] = meshgrid(1:size(cone,2),1:size(cone,1));
         E_coor = E_coor - dirac_Es(jjj);  
         K_coor = K_coor - dirac_ks(jjj);  
         interp_arpes = interp2(E_coor, K_coor, cone, E_interp, K_interp);
         regional_arpes(:,:,iii) = regional_arpes(:,:,iii) + interp_arpes;
         
         DPE_pix = dirac_Es(jjj);
         FL_pix = rfc_FL_Es_(jjj);
         DPE_eV = (-FL_pix + DPE_pix) * pix2eV;
         FL_eV = FL_pix * pix2eV;
         mean_DPE_eV(iii) = mean_DPE_eV(iii) + DPE_eV;
         mean_FL_eV(iii) = mean_FL_eV(iii) + FL_eV;
     end
  
         
     mean_DPE_eV(iii) = mean_DPE_eV(iii) / length(find(region_list{iii}==1));
     mean_FL_eV(iii) = mean_FL_eV(iii) / length(find(region_list{iii}==1));
     nnn_scans(iii) = nnn;
     NNN = NNN + nnn;
end


% plot out the binned and symmetrized regional spectra

E_bin=10;
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
for iii=1:size(region_list,1)    
        
    region_arpes_binned = Binning_2d(regional_arpes(:,:,iii),E_bin,K_bin);
    [region_arpes_symmetrized,Kaxis]=Symmetrized_spectra(region_arpes_binned,K_binned);
    ax = subplot(2,size(B_interval_list,1),iii); %%%%%%replace 2 with size(B_interval_list,1)%%%
    
    %imagesc(regional_arpes(:,:,iii)), axis xy
    %imagesc(rot90(Eaxis,-1),rot90(Kaxis,-1),rot90(region_arpes_binned(1:103,:),-1)), %axis xy;
    %imagesc(Eaxis,Kaxis,region_arpes_symmetrized), axis xy;
    Eaxis_eV = E_interp(1,1:E_bin:end-1) * pix2eV;
    Eaxis_eV_0fixed = Eaxis_eV - abs(mean_DPE_eV(iii));           
    Kaxis_invA = Kaxis * pix2invA;
    if length(find(region_list{iii}==1)) > 0
        %imagesc(Kaxis_invA, flip(Eaxis_eV_0fixed), rot90(norman(region_arpes_symmetrized,0,5))), axis xy
        imagesc(Kaxis_invA, flip(Eaxis_eV_0fixed), rot90(norman(region_arpes_binned(1:103,:),0,5))), axis xy

        colormap hot
        %yticks([0,-(mean_DPE_eV(iii)/pix2eV)])
        %yticklabels({num2str(mean_DPE_eV(iii)),'0'})
        ylabel_str = 'E - E$_{F}$, (eV)';
        ylabel(ylabel_str,'Interpreter','latex','FontSize',8)
        %yt = get(gca, 'YTick');
        %set(gca, 'FontSize',7)
        %ytickangle(90)
        ylim([-0.75,0.1])
        %xticks([-50,0,50])
        %xticklabels({num2str(-50*pix2invA),'0',num2str(50*pix2invA)})
        xlabel('Momentum ({\AA}$^{-1}$)','Interpreter','latex','FontSize',8)
        %xt = get(gca, 'XTick');
        %set(gca, 'FontSize', 7)
    end
    if iii == 1+size(B_interval_list,1)*nn
        descrr = {[num2str(A_interval_list(nn+1,1)),'-',num2str(A_interval_list(nn+1,2))]};
        axes(ax)
        text(-0.2,-0.3,descrr')
        %ylabel([num2str(A_interval_list(nn+1,1)),'-',num2str(A_interval_list(nn+1,2))]);
        nn=nn+1;
    end 
    if ismember(iii,[1:size(B_interval_list,1)]) == 1
        %set(gca,'XaxisLocation','top');
        descr = {[num2str(B_interval_list(nnn+1,1)),'-',num2str(B_interval_list(nnn+1,2))]};
        axes(ax)
        text(-.05,.15,descr)
        %xlabel([num2str(B_interval_list(nnn+1,1)),'-',num2str(B_interval_list(nnn+1,2))]);
        nnn=nnn+1;
    end
    ax2 = subplot(2,size(B_interval_list,1),iii+size(B_interval_list,1));
    imagesc(E_small_map), axis xy, hold on;
    %plot(map_regions(1,:,iii),map_regions(2,:,iii),'b.'), hold off;
    title(['nnn=',num2str(nnn_scans(iii))])
end
suptitle([table_title,', total scans: ',num2str(NNN)]);




