tic;
load_data = 0;
if load_data == 1
    load 'cones.mat'
    load 'kLOS.mat'
    load 'rfc_big_scan_170927.mat';
    load 'cones55555.mat';
    load rfc_FL_scan_170927.mat;
end

%%%%%PARAMS%%%%%%%%%%%%%%%%
scan_is =   [round(961*rand), round(961*rand), round(961*rand), round(961*rand)];%(Neighbor_sites(12,14,31,31))';%round(961*rand);% [round(961*rand), round(961*rand), round(961*rand), round(961*rand)];

%bcb_drawfit;

wannasee = 1;
wannasee_B0Bs = 0;

bin_E = 5;
bin_k = 2;

fass_sigma = 5;
%b_E_1 = 25;
%b_E_2 = 105;
%b_k_1 = 40;
%b_k_2 = 140;

IT_sigma = 2;

LineWidth = 1;
draw_sigma = 4;

padding = 2*5;
draw_box = zeros([81+padding],[41+padding]);
E_0 = 40 + padding/2;
K_0 = round(size(draw_box,1)/2);
draw_x = (1:size(draw_box,1))';
draw_x0 = draw_x - K_0;

B0_range = [5:1:30];%1.0];%[1.3:.05:2.0];%[1.3:.1:2.5];%[0.5:.5:3.0];%[.6:.2:2];%[1.8];%[1,2,3];
B_range =  [.03:.005:.08];%[0:.002:.01];%[0,.05,.1];

%B0_range_eVA = B0_range * (bin_E) * pix2eV;
B_range_eVA = B_range * (bin_E/bin_k) * .8107;

corrspread_TH = 0.95;
%%%%%%%%%%%%%%%%%%%%
%{
Es_ato = zeros(1,num_scans);
ks_ato = zeros(1,num_scans);
As_ato = zeros(1,num_scans);
Bs_ato = zeros(1,num_scans);
corrs_ato = zeros(1,num_scans);
corrspreads_ato = zeros(1,num_scans);
%}
%First draw all the A,B templates%%%
B0B_ITs = zeros(size(draw_box,1)-padding,size(draw_box,2)-padding,length(B0_range)*length(B_range));
%ABBA_IT10s = zeros(size(draw_box,1)-padding,size(draw_box,2)-padding,length(B0_range)*length(B_range));
nan_B0Bs = []; ABBA_IT_i = 1; nan_AB_i = 1;
if wannasee_B0Bs == 1
    figure
end

for B0_i = 1:length(B0_range)
    B0 = B0_range(B0_i);
    for B_i = 1:length(B_range)
        B = B_range(B_i);

        draw_yp = B*((draw_x0).^2) + B0;

        draw_itp = horzcat(draw_yp,draw_x);

        %draw_itp_curt = draw_itp(max(1,round(-B0/(2*B)+K_0)):end,:);
        draw_itp_curt = draw_itp;
        draw_itp_it = reshape(draw_itp_curt',1,[]);
        
        if length(B_range) == 1
            draw_itp_it_ = zeros(1,4);
            draw_itp_it_(1:2) = draw_itp_it(1:2);
            draw_itp_it_(3:4) = draw_itp_it(end-1:end);
            draw_itp_it = draw_itp_it_;
        end      
        
        ITP = insertShape(draw_box, 'Line', draw_itp_it, 'LineWidth',LineWidth);

        ITP = ITP(:,:,1)+ITP(:,:,2) + ITP(:,:,3);
        
        ITT = ITP;
        ITT(ITT~=0) = 1;
        
        ITT_ = abs(1 - ITT);
        %figure, imagesc(ITT_), axis xy, title('ITT\_')
        wedge_labels = bwlabel(ITT_);
        label_bot = wedge_labels(round(size(ITT_,1)/2),1);
        label_top = wedge_labels(round(size(ITT_,1)/2),end);
        %label_left = wedge_labels(1,round(size(ITT_,2)/2));
        %label_right = wedge_labels(end,round(size(ITT_,2)/2));
        
        if label_top ~= 0
            ITT(wedge_labels == label_top) = 1;
        end
        if label_bot ~= 0
            ITT(wedge_labels == label_bot) = 0;
        end
%         if label_left ~= 0
%             ITT(wedge_labels == label_left) = NaN;
%         end
%         if label_right ~= 0
%             ITT(wedge_labels == label_right) = NaN;
%         end
        
        %figure, imagesc(ITT), axis xy, title('ITT')
        
        IT = imgaussfilt( ITT, IT_sigma);
        %figure, imagesc(IT), axis xy, title('IT')
        ITtt = zeros(size(IT));
        for ititi = 1:size(IT,1)
            ITtt(ititi,:) = mat2gray(IT(ititi,:));
        end
        %ITcgm = ITcg - mean(ITcg(:));
        ITc = ITtt((padding/2)+1:end-(padding/2),(padding/2)+1:end-(padding/2));
        ITcg = mat2gray(ITc);
        ITcg(ITcg==0) = NaN;
%         IT_ = ITcg;
%         IT_((ITc~=0))=1;
%         IT_ = abs(1-IT_);
        
        
        
        %bot_top_wedge_finder = ITcg(round(size(ITcg,1)/2),:);
        
        
        ITcgc = ITcg;%IT_processor(ITcg,0,0);
        ITcgm = ITcg;%c - nanmean(ITcgc(:));
        %IT10 = ITcgm;
        %IT10(isnan(IT10)==0) = 1;
        %IT10(isnan(IT10)) = 0;
        
        if wannasee_B0Bs == 1       
            
            AAA_i = ceil(ABBA_IT_i/length(B_range));
            BBB_i = ABBA_IT_i - (AAA_i-1)*length(B_range);

            show_IT = ITcgm;%((padding/2)+1:end-(padding/2),(padding/2)+1:end-(padding/2));
            show_IT(isnan(show_IT)) = 100;
            subplot(length(B0_range),length(B_range),...
            (length(B0_range)-AAA_i)*length(B_range) + BBB_i)
            imagesc(show_IT), axis xy, axis off;
            pause(.001)
        end
        
        B0B_ITs(:,:,ABBA_IT_i) = ITcgm;
        %ABBA_IT10s(:,:,ABBA_IT_i) = IT10;
        ABBA_IT_i = ABBA_IT_i+1;
    end
end
disp('Finished drawing templates; '), toc

%%%For each Scan, do ssd scan with each of the ABBA_ITs(:,:,i)%%%
for i = scan_is
    if rem(i,50) == 0
        disp(['Starting scan #',num2str(i)]), toc
    end
    scone = cones(:,:,i);
    scone_f = imgaussfilt(scone,fass_sigma);   
    scone_bf = Binning_2d(scone_f, bin_E, bin_k);  
    b_k_1 = round(kLOS(i)/bin_k)-40;
    b_k_2 = round(kLOS(i)/bin_k)+40;
    b_E_1 =  round(rfc_FL_Es(i)/bin_E + bin_E/2)-50;
    b_E_2 = round(rfc_FL_Es(i)/bin_E + bin_E/2);
    scone_bfc = scone_bf(b_k_1:b_k_2, b_E_1:b_E_2);
    fass = scone_bfc; 
    fass_g = mat2gray(fass);
    %fass_gm = fass_g - mean(fass_g(:));
    fass_gg = fass_g;
    fass_gg(fass_gg>(mean(fass_g(:))+.5*std(fass_g(:)))) = mean(fass_g(:))+.5*std(fass_g(:));
    f_image = fass_gg; 
    
    %figure, imagesc(scone_bfc), axis xy
    
    tempsize = size(B0B_ITs(:,:,1));
    B0B_corr_table = zeros([length(B0_range),length(B_range)]);  
    corr_IT = zeros([(size(f_image,1)-tempsize(1)+1),(size(f_image,2)-tempsize(2)+1)]);

    for AB_IT_i = 1:(length(B0_range)*length(B_range))
        B0_i = ceil(AB_IT_i / length(B_range));
        B_i = AB_IT_i - (B0_i-1)*length(B_range);
        if ismember(AB_IT_i, nan_B0Bs) == 1      
            B0B_corr_table(B0_i,B_i) = NaN;
            continue
        end
        
        IT_template = B0B_ITs(:,:,AB_IT_i);
        IT_mid = IT_template(round(size(IT_template,1)/2),:);
        IT_weight = 1;%/(1-0.002*(size(IT_template,2) - find(IT_mid == 1,1, 'first')));
        
        %IT10_temp = ABBA_IT10s(:,:,AB_IT_i);
        N = size(IT_template,1);
        M = size(IT_template,2);

        corr = zeros(size(corr_IT));
        
        for x = 1%:(size(f_image,2) - size(IT_template,2) +1) 
            for y = 1:(size(f_image,1) - size(IT_template,1) +1)  
                f_image_under = f_image(y:y+N-1, x:x+M-1);
                corr(y,x) =  IT_weight*sum(nansum( IT_template .* f_image_under ))/nansum(nansum(IT_template));
            end
        end
        [corrmax_row,corrmax_col] = find(corr==max(corr(:)));
        B0B_corr_table(B0_i,B_i) = max(corr(:));
        
        if B0B_corr_table(B0_i,B_i) == max(B0B_corr_table(:))
            disp('updating its')
            corr_IT = corr;
            B0_i_it = B0_i;
            B_i_it = B_i;
            AB_IT_i_it = AB_IT_i;
        end 
    end
    %figure, imagesc(ABBA_ITs(:,:,AB_IT_i_it)), axis xy
    [corrmax_row_it,corrmax_col_it] = find(corr_IT==max(corr_IT(:)));
    corr_ITg = (corr_IT - min(corr_IT(:)))/(max(corr_IT(:))-min(corr_IT(:)));
    corrspread = length(find(corr_ITg > corrspread_TH));

    DPE_b = B0_range(B0_i_it) + b_E_1 - 1;
    DPK_b = K_0 + b_k_1 - 1;
    DPE_pix = (DPE_b * bin_E) + bin_E/2;
    DPK_pix = (DPK_b* bin_k) + bin_k/2;

    Es_ato(i) = DPE_pix;
    ks_ato(i) = DPK_pix;
    B0s_ato(i) = B0_range(B0_i_it);
    Bs_ato(i) = B_range(B_i_it);
    corrs_ato(i) = max(corr_IT(:));
    corrspreads_ato(i) = corrspread;

    %figure, subplot(1,2,1), imagesc(f_image), axis xy, subplot(1,2,2), imagesc(overlay), axis xy


    if wannasee == 1
        %overlay = zeros(size(f_image));
        overlay = zeros(size(scone_bf));
        overlay((1:size(IT_template,1))+b_k_1-1,(1:size(IT_template,2))+b_E_1-1) = B0B_ITs(:,:,AB_IT_i_it);
        overlay(isnan(overlay)) = -.1;
        
        figure, 
        subplot(2,2,1), 
        imagesc(B_range,B0_range,B0B_corr_table), axis xy, hold on;
        plot(B_range(B_i_it),B0_range(B0_i_it),'w+'), hold off
        title({['mc=',num2str(corrs_ato(i))];['A=',num2str(B0_range(B0_i_it)),' B=',num2str(B_range(B_i_it))]})

        subplot(2,2,2),
        imagesc((1:size(corr_IT,1)),fliplr(1:size(corr_IT,2)),rot90(corr_IT,1)), axis xy, hold on
        plot(corrmax_row_it,corrmax_col_it,'w+'), hold off
        title({['spread=',num2str(corrspreads_ato(i))];['maxC=',num2str(corrs_ato(i))]})

        subplot(2,2,3), 
        %imagesc((1:size(f_image,1)),fliplr(1:size(f_image,2)),rot90(f_image,1)), axis xy
        imagesc(scone_bf'), axis xy, hold on;
        plot([0,size(scone_bf,2)],[DPE_b,DPE_b],'r'), hold off;
        title(['DPI=',num2str(round(DPI_big(i)))])
        
        
        %subplot(2,2,4), imagesc((1:size(f_image,1)),fliplr(1:size(f_image,2)),rot90(f_image+.5*overlay,1)), axis xy, hold on;
        %plot([1,size(IT_template,1),size(IT_template,1),1,1]+corrmax_row_it-1,[1,1,size(IT_template,2),size(IT_template,2),1]+corrmax_col_it-1,'w'), hold off;
        subplot(2,2,4), imagesc((1:size(scone_bf,1)),fliplr(1:size(scone_bf,2)),rot90(scone_bf+5*overlay,1)), axis xy, hold on;
        plot([1,size(IT_template,1),size(IT_template,1),1,1]+corrmax_row_it-1+b_k_1-1,[1,1,size(IT_template,2),size(IT_template,2),1]+corrmax_col_it-1+b_E_1-1,'w'), hold off;
        title(['Epix=',num2str(DPE_pix),' Kpix=',num2str(DPK_pix)])

        suptitle(['i=',num2str(i)])
        
        
        pause(.001)
    end
    
end

toc;

%}


