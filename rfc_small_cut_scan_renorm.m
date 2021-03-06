tic; 

rfc_small_cut_Es = zeros(1,num_scans);
rfc_small_cut_ks = zeros(1,num_scans);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bin_E = 3;
bin_k = 2;

rfc_cut_sigma = 3;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Making the 0th-Order Template %%%%%%%%%%%%%%%%%%%%%%%%%%%
%{
raw_full_cone = zeros(size(cone_range_K,2),size(cone_range_E,2));

for i = 1:1:num_scans;
    raw_full_cone = raw_full_cone + cones(:,:,i);
end
%}
rfc_b = Binning_2d(raw_full_cone, bin_E, bin_k);
rfc_b_zeros = zeros(size(rfc_b));
rfc_b_cutter = insertShape(rfc_b_zeros,'Line',[195,72, 146,90, 195,106]);
rfc_b_cutter = rfc_b_cutter(:,:,1);
rfc_b_cutter(rfc_b_cutter~=0) =1;
for r_i = 146:195
    rfc_b_col = rfc_b_cutter(:,r_i);
    pntA = find(rfc_b_col==1,1,'first');
    pntB = find(rfc_b_col==1,1,'last');
    rfc_b_col(pntA:pntB) = 1;
    rfc_b_cutter(:,r_i) = rfc_b_col;
end
rfc_b_cutter = abs(rfc_b_cutter-1);
rfc_b_cut = rfc_b_cutter.*rfc_b;
rfc_small_cut_window = rfc_b_cut(50:130,80:160);
rfc_small_cut_bwg = mat2gray(rfc_small_cut_window);

window_small_cut_k_coor = 42;
window_small_cut_E_coor = 28; %%%%%%%%%%%%%%% hand-pick point from rfc_bg%%%%%%%%%%%

wedge = rfc_small_cut_bwg;
wedge(wedge==0) = 100;
wedge(wedge~=100) = 0;
wedge = wedge ./ 100; 
%wedge(:,end-5:end) = 0;

disp('Template step completed '), toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Round 1 Scan with Large Window Including Bulk Bands %%%%%%%%%%%%%

scan_maxes_small_cut = zeros(1,num_scans);
%DPI_small_cut = zeros(1,num_scans);
%scan_spread_cut = zeros(1,num_scans);
wedge_sum = zeros(1,num_scans); 

NNN=0;
for i = [round(961*rand),round(961*rand),round(961*rand),round(961*rand),round(961*rand), round(961*rand)]%1:1:num_scans;  
    NNN=NNN+1;
    if rem(i,100) == 0
        disp(['Now on scan ',num2str(i)])
    end
    
    scone = result(:,:,i);
    Epix = (rfc_big_Es(i) - bin_E/2)/bin_E;
    kpix = (rfc_big_ks(i) - bin_k/2)/bin_k;
    Epixhr = size(rfc_small_cut_bwg,2)/2 + 10; 
    kpixhr = size(rfc_small_cut_bwg,1)/2 + 2;
    
    scone = imgaussfilt(scone, rfc_cut_sigma);
    scone_b = Binning_2d(scone, bin_E, bin_k);
    scone_bc = scone_b( max(1,round(kpix-kpixhr)):min(end,round(kpix+kpixhr)), max(1,round(Epix-Epixhr)):min(end,round(Epix+Epixhr+5)));                    
    scone_bg = mat2gray(scone_bc);
    
    rfc_corr_table = zeros((size(scone_bg,1)-size(rfc_small_cut_bwg,1)+1), (size(scone_bg,2)-size(rfc_small_cut_bwg,2)));
    E_i = 0; 
    for E_off = (0 : (size(scone_bg,2) - size(rfc_small_cut_bwg,2)))
        E_i = E_i +1;
        K_i = 0;
        for K_off = (0 : (size(scone_bg,1) - size(rfc_small_cut_bwg,1)))
            K_i = K_i + 1;
            scone_window = scone_bg((1:size(rfc_small_cut_bwg,1)) + K_off, (1:size(rfc_small_cut_bwg,2)) + E_off);
            rfc_corr_table(K_i, E_i) = sum(dot(rfc_small_cut_bwg, scone_window)) / sum(dot(abs(wedge-1),scone_window));
        end
    end
    
    rfc_corr_table_norm = (rfc_corr_table - min(rfc_corr_table(:)))/(max(rfc_corr_table(:))-min(rfc_corr_table(:)));
    rfc_corr_spread = length(find(rfc_corr_table_norm > 0.9))/(size(rfc_corr_table_norm,1)*size(rfc_corr_table_norm,2));
    
    [rfc_corr_max, rfc_corr_max_i ] = max(rfc_corr_table(:));
    [rfc_K_off, rfc_E_off] = ind2sub(size(rfc_corr_table), rfc_corr_max_i);
    
    scone_bc_DPK = rfc_K_off + window_small_cut_k_coor;
    scone_bc_DPE = rfc_E_off + window_small_cut_E_coor;
    
    rfc_DPK_b = scone_bc_DPK + max(1,round(kpix-kpixhr)) -1;
    rfc_DPE_b = scone_bc_DPE + max(1,round(Epix-Epixhr)) -1;
    
    rfc_DPK = (rfc_DPK_b)*bin_k + (bin_k/2);
    rfc_DPE = (rfc_DPE_b)*bin_E + (bin_E/2);
    
    overlay = zeros(size(scone_b));
    overlay((1:size(rfc_small_cut_bwg,1))+rfc_DPK_b - window_small_cut_k_coor,(1:size(rfc_small_cut_bwg,2))+rfc_DPE_b - window_small_cut_E_coor) = rfc_small_cut_bwg;
    overlay = mat2gray(overlay);
    overlay(overlay < 0.3) = 0;
    wedge_overlay = zeros(size(scone_b));
    wedge_overlay((1:size(rfc_small_cut_bwg,1))+rfc_DPK_b - window_small_cut_k_coor,(1:size(rfc_small_cut_bwg,2))+rfc_DPE_b - window_small_cut_E_coor) = wedge;
    
    wedge_sum(i) = sum(dot(wedge_overlay, mat2gray(scone_b)));
    
    rfc_small_cut_ks(i) = rfc_DPK;
    rfc_small_cut_Es(i) = rfc_DPE;
    scan_maxes_small_cut(i) = rfc_corr_max;
    
    
    
    %{
    figure
    imagesc(mat2gray(scone_b)+.3*wedge_overlay), axis xy, hold on;
    plot(rfc_DPE_b, rfc_DPK_b, 'r*'), hold off;
    %}
    ix = imregionalmax(rfc_corr_table);
    [ixr,ixc] = find(ix==1);
    
    figure
    subplot(2,2,[1 3])
    %imagesc(rot90(rfc_corr_table,-1)), axis xy;
    contourf(rot90(rfc_corr_table,-1),50), hold on;
    plot(size(rfc_corr_table,1)-ixr+1,ixc,'r*');
    title({['MaxC=',num2str(round(rfc_corr_max,3))];['spread=',num2str(round(rfc_corr_spread,3))]},'FontSize',8)

    subplot(2,2,2)
    E = mat2gray(scone_b);
    Ea = rot90(E,-1);
    I = mat2gray(overlay);
    Ia = rot90(I,-1);
    greenlol = cat(3,.7*ones(size(Ea)),zeros(size(Ea)),zeros(size(Ea)));
    imshow(Ea,'InitialMag','fit','DisplayRange',[0 .5], 'XData',[0 1], 'YData', [0 1]), hold on;
    h = imshow(greenlol, 'XData',[0,1], 'YData',[0,1]);  axis xy
    set(h,'AlphaData',Ia);
    hold off;
    title({['i=',num2str(i)]; ['ws=',num2str(wedge_sum(i))]})
    
    subplot(2,2,4)
    imagesc(rot90(scone_b,-1)), axis xy
    title({['MC*DPI=',num2str(round(rfc_corr_max*DPI_big(i),3))];...
        ['E=',num2str(rfc_DPE),' K=',num2str(rfc_DPK)]})
    
    %{
    scann = normxcorr2( rfc_small_cut_bwg, scone_bg );
    scan=scann(size(rfc_small_cut_bwg,1)+1:end-size(rfc_small_cut_bwg,1),size(rfc_small_cut_bwg,2)+1:end-size(rfc_small_cut_bwg,2));
    [scan_max, scan_max_i] = max(scan(:));
    [smax_y,smax_x] = ind2sub(size(scan),scan_max_i);
    scan_maxes_small_cut(i) = scan_max;
    
    [y_peak, x_peak] = find(scan==max(scan(:)));
    y_offset = y_peak - size(rfc_small_cut_bwg,1) + size(rfc_small_cut_bwg,1);
    x_offset = x_peak - size(rfc_small_cut_bwg,2) + size(rfc_small_cut_bwg,2);
        
    rfc_small_cut_Es(i) = (x_offset + window_small_cut_E_coor + round(Epix-Epixhr) -1)*bin_E + bin_E/2; %%%%%%!!!!%%
    rfc_small_cut_ks(i) = (y_offset + window_small_cut_k_coor + round(kpix-kpixhr) -1)*bin_k + bin_k/2;    
   
    scan_gray = mat2gray(scan);
    
    wedge_sum(i) = sum(dot(wedge(:,1:end-5), scone_bg(y_offset+1:y_offset+size(rfc_small_cut_bwg,1),x_offset+1:x_offset+size(rfc_small_cut_bwg,2)-5)));
    
    
    if (smax_y-5 <= 1) || (smax_y+5>=size(scan_gray,1)) || (smax_x-5 <= 1) || (smax_x+5 >= size(scan_gray,2))
        scan_spread_cut(i) = NaN;
        DPI_small_cut(i) = NaN;
        %continue
    else
        %scan_spread_cut(i) = sum(sum(scan_gray(smax_y-5:smax_y+5,smax_x-5:smax_x+5)));
        scan_spread_cut(i) = sum(scan(:));
        DPI_small_cut(i) = sum(sum(scone_bg(max(1,y_offset):min(size(scone_bg,1),y_offset+size(rfc_small_cut_bwg,1)),window_small_cut_E_coor+x_offset-10:window_small_cut_E_coor+x_offset+10)));
    end
    
    
    figure;
    subplot(1,4,1), %imagesc(scan_Kaxis,scan_Eaxis,rot90(scan,-1)), axis xy, title(num2str(i)), hold on;
    imagesc(rot90(scan,-1)), axis xy, hold on;
    plot(size(scan,1)-y_peak+1,x_peak,'b*'), hold off;
        title({['Sc max=',num2str(round(scan_maxes_small_cut(i),3))];['spread=',num2str(round(scan_spread_cut(i)))]});
        
    subplot(1,4,2), imagesc(rot90(rfc_small_cut_bwg,-1)), axis xy, hold on;
    plot([window_small_cut_k_coor],[window_small_cut_E_coor],'r*'), hold off;
    title(['ws=',num2str(round(wedge_sum(i),1))]);

    
    subplot(1,4,3), imagesc(rot90(scone_bg,-1)), axis xy, hold on;
    %plot([1,size(rfc_small_bwg,2),size(rfc_small_bwg,2),1,1]+x_offset,[1,1,size(rfc_small_bwg,1),size(rfc_small_bwg,1),1]+y_offset,'r'), hold off;
    plot([0,0,-size(rfc_small_cut_bwg,1),-size(rfc_small_cut_bwg,1),1]+size(scone_bg,1)-y_offset,[1,size(rfc_small_cut_bwg,2),size(rfc_small_cut_bwg,2),1,1]+x_offset,'r'), hold on;
    plot([size(scone_bg,1)-y_offset-window_small_cut_k_coor],[x_offset+window_small_cut_E_coor],'r*'), hold off;
    %plot(flipud([1,size(rfc_small_bwg,1),size(rfc_small_bwg,1),1,1]+size(scone_bg,1)-size(rfc_bwg,1)+y_offset),[1,1,size(rfc_bwg,2),size(rfc_bwg,2),1]+x_offset,'r'), hold off;
    title(['E=',num2str(rfc_small_cut_Es(i)),' k=',num2str(rfc_small_cut_ks(i))]);
    
    subplot(1,4,4), imagesc(rot90(imgaussfilt(Binning_2d(cones(:,:,i),3,2),rfc_cut_sigma),-1)), axis xy;
    title(['i=',num2str(i)]);
    %}
    
end
disp('Scanning Round 1 completed'), toc


E_small_map = reshape(rfc_small_cut_Es, 31, 31);
k_small_map = reshape(rfc_small_cut_ks, 31, 31);

EE_small_sd = mean([std(E_small_map(:,2:end)-E_small_map(:,1:end-1)),std(E_small_map(2:end,:)-E_small_map(1:end-1,:))])
kk_small_sd = mean([std(k_small_map(:,2:end)-k_small_map(:,1:end-1)),std(k_small_map(2:end,:)-k_small_map(1:end-1,:))])

%{
figure, subplot(2,2,[1 2]), imagesc(rfc_bwg), axis xy, title('window1');
subplot(2,2,3), imagesc(E_map), axis xy, title('Emap1');
subplot(2,2,4), imagesc(k_map), axis xy, title('kmap1');

figure, subplot(1,2,1), imagesc(E_map), axis xy, title('Emap');
subplot(1,2,2), imagesc(k_map), axis xy, title('kmap');
%}
