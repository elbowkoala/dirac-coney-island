tic; 
wannasee = 0;
scan_is = 1:961;

rfc_Es_before = zeros(1,num_scans);
rfc_ks_before = zeros(1,num_scans);
rfc_wedges_before = zeros(1,num_scans); 
rfc_corrspreads_before = zeros(1,num_scans);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin_E = 3;
bin_k = 2;

ssd_sigma = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Making the 0th-Order Template %%%%%%%%%%%%%%%%%%%%%%%%%%%

raw_full_cone = zeros(size(cone_range_K,2),size(cone_range_E,2));

for i = 1:1:num_scans;
    raw_full_cone = raw_full_cone + cones(:,:,i);
end


rfc_b = Binning_2d(raw_full_cone, bin_E, bin_k);
rfc_patch = rfc_b(51:130,81:160);
%rfc_patch = rfc_b(40:140,65:210);
%rfc_patch = rfc_b(:,195:235);
%[rfc_patch_sym, rfc_patch_symkaxis] = Symmetrized_spectra(rfc_patch,(1:size(rfc_patch,1))');
rfc_patch_gray = mat2gray(rfc_patch);
ssd_template = rfc_patch_gray - mean(rfc_patch_gray(:));

template_DPK = 42;
template_DPE = 37;
%template_DPK = 52;
%template_DPE = 19;

temp_W = size(ssd_template,2);
temp_H = size(ssd_template,1);

wedge = insertShape(zeros(size(ssd_template)),'Line',[80,32, 60,40, 80,48]);
wedge = wedge(:,:,1);
wedge(wedge~=0) = 1;
for w_i = 60:size(ssd_template,2)
    wedge_col = wedge(:,w_i);
    pntA = find(wedge_col==1,1,'first');
    pntB = find(wedge_col==1,1,'last');
    wedge_col(pntA:pntB) = 1;
    wedge(:,w_i) = wedge_col;
end

mean_ssd_big_k = mean(ssd_big_ks);

tic
for i = scan_is%[round(961*rand),round(961*rand),round(961*rand),round(961*rand),round(961*rand), round(961*rand)]%1:1:num_scans;  
    if rem(i,100) == 0
        disp(['Now on scan ',num2str(i), ';  time: ',num2str(toc)])
    end
    b_E = round((ssd_big_Es(i)-bin_E/2)/bin_E);
    b_k = round((mean_ssd_big_k-bin_k/2)/bin_k);
    scone_b_E_range = [round( b_E - size(ssd_template,2)/2 - 25) : round(b_E + size(ssd_template,2)/2 + 35)];
    scone_b_k_range = [round( b_k - size(ssd_template,1)/2 - 15) : round(b_k + size(ssd_template,1)/2 + 15)];
    
    
    
    scone = cones(:,:,i);
    
    scone = imgaussfilt(scone, ssd_sigma);
    scone_b = Binning_2d(scone, bin_E, bin_k);
    scone_b = scone_b(scone_b_k_range,scone_b_E_range);
    %scone_b = scone_b(:,171:260);
    
    scone_bg = mat2gray(scone_b);
    scone_bgm = scone_bg - mean(scone_bg(:));
    
    ssd_image = scone_bgm;
    
    corr = zeros([(size(ssd_image,1)-size(ssd_template,1)+1),(size(ssd_image,2)-size(ssd_template,2)+1)]);
    
    for x = 1:(size(ssd_image,2) - size(ssd_template,2) +1)
        for y = 1:(size(ssd_image,1) - size(ssd_template,1) +1)
            for m = 1:size(ssd_template,2)
                for n = 1:size(ssd_template,1)
                    corr(y,x) = corr(y,x) + (ssd_image(y+n-1, x+m-1) - ssd_template(n,m)).^2;
                end
            end
        end
    end
    
    [ssdmin_row,ssdmin_col] = find(corr==min(corr(:)));
    

    ssd_E_b = ssdmin_col + template_DPE + scone_b_E_range(1)-1 - 1;
    ssd_k_b = ssdmin_row + template_DPK +scone_b_k_range(1)-1 - 1;

    rfc_Es_before(i) = ssd_E_b * bin_E + bin_E/2;
    rfc_ks_before(i) = ssd_k_b * bin_k + bin_k/2;
    
    
    overlay = zeros(size(ssd_image));
    wedge_overlay = zeros(size(ssd_image));
    overlay((1:size(ssd_template,1))+ssdmin_row-1,(1:size(ssd_template,2))+ssdmin_col-1) = ssd_template;
    wedge_overlay((1:size(ssd_template,1))+ssdmin_row-1,(1:size(ssd_template,2))+ssdmin_col-1) = wedge;
    rfc_wedges_before(i) = sum(dot(wedge_overlay,scone_bg))/(nnz(wedge_overlay));
    [wedge_row,wedge_col] = find(wedge~=0);
    
    corr_gray = mat2gray(corr);
    rfc_corrspreads_before(i) = length(find(corr_gray<0.1))/(size(corr,1)*size(corr,2));
    
    if wannasee == 1
        
        figure

        subplot(2,2,1)
        imagesc((1:size(ssd_template,1)), fliplr(1:size(ssd_template,2)),rot90(ssd_template,1)), axis xy, hold on;
        plot([wedge_row],[wedge_col],'r.'), hold off
        title(['mincorr = ',num2str(round(min(corr(:))))])


        subplot(2,2,2)
        imagesc((1:size(ssd_image,1)),fliplr(1:size(ssd_image,2)),rot90(ssd_image,1)), axis xy
        title(['DPI = ',num2str(round(DPI_big(i)))])

        subplot(2,2,3)
        imagesc((1:size(corr,1)),fliplr(1:size(corr,2)),rot90(corr,1)), axis xy, hold on
        plot(ssdmin_row,ssdmin_col,'w+'), hold off
        title(['spread = ',num2str(round(rfc_corrspreads_before(i),3))])


        subplot(2,2,4)
        imagesc((1:size(ssd_image,1)),fliplr(1:size(ssd_image,2)),rot90(overlay+ssd_image,1)), axis xy
        title(['ws = ',num2str(round(rfc_wedges_before(i),3))])

        suptitle(['i=',num2str(i),'  E=',num2str(rfc_Es_before(i)),'  k=',num2str(rfc_ks_before(i))])
        pause(.01)
    end
end
%disp('Scanning Round 1 completed'), toc
toc

E_small_map = reshape(rfc_Es_before, 31, 31);
k_small_map = reshape(rfc_ks_before, 31, 31);

EE_small_sd = mean([std(E_small_map(:,2:end)-E_small_map(:,1:end-1)),std(E_small_map(2:end,:)-E_small_map(1:end-1,:))]);
kk_small_sd = mean([std(k_small_map(:,2:end)-k_small_map(:,1:end-1)),std(k_small_map(2:end,:)-k_small_map(1:end-1,:))]);

