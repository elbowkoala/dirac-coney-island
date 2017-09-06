

table_title = 'scan 170829b, Es vs ks';
A_map = reshape(draw_Es_905_tails,31,31);      %Vertical axis 
B_map = reshape(draw_ks_905_tails,31,31);%reshape(line_ks_829,31,31); %reshape(ABEK_ks,31,31);%DPI_map;   %Horizontal axis
dirac_Es = draw_Es_905_tails;   %input Energies vector (cone pixels)
dirac_ks = draw_ks_905_tails;   %input mtm vector (cone pixels)

A_interval_list=[0,1; .3,.7; .4,.65];%; .3,.7; .2,.5; .5,.8];%[0,1; 0,0.3; 0.6,1; .2,.8];% 0,0.2; .1,.9; .3,.7]; %Setting the intervals 
B_interval_list = [0,1; .3,.7; .4,.65];%[0,1; 0,0.3; 0.6,1; .2,.8];% 0,0.8; .2,.8; 0,.7]; %[0,1; .3,1; .4,1; .5,1];%[0,1; .1,.9; .3,.7];%[0.1,1; 0.3,1; 0.6,1]; 




[E_interp,K_interp] = meshgrid(-149.5:.5:300, -80:.5:80); %meshgrid(-149.5:.5:280, -80:.5:80); %meshgrid(-199.5:0.5:250,-80:0.5:80);
full_arpes = zeros(size(E_interp));


[region_list]=Correlation_list_indep_ranges(A_map,B_map,A_interval_list,B_interval_list);
regional_arpes=zeros([size(E_interp),size(region_list,1)]);
E_coor_ave = zeros(300,800,9);
K_coor_ave = zeros(300,800,9);

for iii=1:size(region_list,1)
    
     for jjj=find(region_list{iii}==1)'
         
         cone = cones(:,:,jjj);
         [E_coor, K_coor] = meshgrid(1:size(cone,2),1:size(cone,1));
         E_coor = E_coor - dirac_Es(jjj);  %%%(bin_E=2)
         K_coor = K_coor - dirac_ks(jjj);  %%%(bin_k=1)
         E_coor_ave(:,:,iii) = E_coor_ave(:,:,iii) + E_coor;
         K_coor_ave(:,:,iii) = K_coor_ave(:,:,iii) + K_coor;
         interp_arpes = interp2(E_coor, K_coor, cone, E_interp, K_interp);
         regional_arpes(:,:,iii) = regional_arpes(:,:,iii) + interp_arpes;
         
     end
     tot_coors = length(find(region_list{iii}==1));
     E_coor_ave(:,:,iii) = E_coor_ave(:,:,iii) / tot_coors;
     K_coor_ave(:,:,iii) = K_coor_ave(:,:,iii) / tot_coors;
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
figure;
nn=0;
nnn=0;
for iii=1:size(region_list,1)    
    region_arpes_binned = Binning_2d(regional_arpes(:,:,iii),E_bin,K_bin);
    [region_arpes_symmetrized,Kaxis]=Symmetrized_spectra(region_arpes_binned,K_binned);
    subplot(size(A_interval_list,1),size(B_interval_list,1),iii)
    %imagesc(regional_arpes(:,:,iii)), axis xy
    imagesc(rot90(Eaxis,-1),rot90(Kaxis,-1),rot90(region_arpes_binned(1:103,:),-1)), %axis xy;
    %imagesc(Eaxis,Kaxis,region_arpes_symmetrized), axis xy;
    if iii == 1+size(B_interval_list,1)*nn
        ylabel([num2str(A_interval_list(nn+1,1)),'-',num2str(A_interval_list(nn+1,2))]);
        nn=nn+1;
    end 
    if ismember(iii,[1:size(B_interval_list,1)]) == 1
        set(gca,'XaxisLocation','top');
        xlabel([num2str(B_interval_list(nnn+1,1)),'-',num2str(B_interval_list(nnn+1,2))]);
        nnn=nnn+1;
    end
end
suptitle(table_title);




