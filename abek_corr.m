tic;
load 'rfc_big_scan_170927.mat';
load 'cones55555.mat';
%%%%%PARAMS%%%%%%%%%%%%%%%%
scan_is =  69%round(961*rand);% [round(961*rand), round(961*rand), round(961*rand), round(961*rand)];%(Neighbor_sites(12,14,31,31))';%round(961*rand);% [round(961*rand), round(961*rand), round(961*rand), round(961*rand)];
wannasee = 1;
wannasee_ABBAs = 0;

bin_E = 5;
bin_k = 2;

fass_sigma = 5;
b_E_1 = 25;
b_E_2 = 105;
b_k_1 = 40;
b_k_2 = 140;

IT_sigma = 2;

LineWidth = 2;
draw_sigma = 4;

padding = 2*5;
draw_box = zeros([71+padding],[61+padding]);
E_0 = 20 + padding/2;
K_0 = round(size(draw_box,1)/2);
draw_x = (1:size(draw_box,1))';
draw_x0 = draw_x - K_0;

A_range = [1.0:.1:2.5];%[1.3:.05:2.0];%[1.3:.1:2.5];%[0.5:.5:3.0];%[.6:.2:2];%[1.8];%[1,2,3];
B_range =  [0:.002:.02];%[0:.002:.01];%[0,.05,.1];

A_range_eVA = A_range * (bin_E/bin_k) * .8107;
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
ABBA_ITs = zeros(size(draw_box,1)-padding,size(draw_box,2)-padding,length(A_range)*length(B_range));
ABBA_IT10s = zeros(size(draw_box,1)-padding,size(draw_box,2)-padding,length(A_range)*length(B_range));
nan_ABs = []; ABBA_IT_i = 1; nan_AB_i = 1;
if wannasee_ABBAs == 1
    figure
end
for A_i = 1:length(A_range)
    A = A_range(A_i);
    for B_i = 1:length(B_range)
        B = B_range(B_i);

        draw_yp = A*(draw_x0) + B*((draw_x0).^2) + E_0;
        draw_yn = -A*(draw_x0) + B*((draw_x0).^2) + E_0;

        draw_itp = horzcat(draw_yp,draw_x);
        draw_itn = horzcat(draw_yn,draw_x);

        draw_itp_curt = draw_itp(max(1,round(-A/(2*B)+K_0)):end,:);
        draw_itn_curt = draw_itn(1:min(length(draw_x),round(K_0+A/(2*B))),:);
        draw_itp_curt = draw_itp;
        draw_itn_curt = draw_itn;
        
        draw_itp_it = reshape(draw_itp_curt',1,[]);
        draw_itn_it = reshape(draw_itn_curt',1,[]);
        
        if length(B_range) == 1
            draw_itp_it_ = zeros(1,4);
            draw_itn_it_ = zeros(1,4);
            draw_itp_it_(1:2) = draw_itp_it(1:2);
            draw_itn_it_(1:2) = draw_itn_it(1:2);
            draw_itp_it_(3:4) = draw_itp_it(end-1:end);
            draw_itn_it_(3:4) = draw_itn_it(end-1:end);
            draw_itp_it = draw_itp_it_;
            draw_itn_it = draw_itn_it_;
        end      
        
        ITP = insertShape(draw_box, 'Line', draw_itp_it, 'LineWidth',LineWidth);
        if nnz(ITP(:,1))==0   %Only consider A,B values that give lines that extend across draw_box
            ABBA_ITs(:,:,ABBA_IT_i) = ABBA_ITs(:,:,ABBA_IT_i) + NaN;
            nan_ABs(nan_AB_i) = ABBA_IT_i;

            ABBA_IT_i = ABBA_IT_i + 1;
            nan_AB_i = nan_AB_i + 1;
            continue
        end
        ITN = insertShape(draw_box, 'Line', draw_itn_it, 'LineWidth',LineWidth);

        ITP = ITP(:,:,1);
        ITN = ITN(:,:,1);
        
        ITT = ITP+ITN;
        ITT(ITT~=0) = 1;
        IT = imgaussfilt( ITT, IT_sigma);
        ITc = IT((padding/2)+1:end-(padding/2),(padding/2)+1:end-(padding/2));
        ITcg = mat2gray(ITc);

        %ITcgm = ITcg - mean(ITcg(:));
        
        IT_ = ITcg;
        IT_((ITc~=0))=1;
        IT_ = abs(1-IT_);
        
        wedge_labels = bwlabel(IT_);
        label_bot = wedge_labels(round(size(IT_,1)/2),1);
        label_top = wedge_labels(round(size(IT_,1)/2),end);
        label_left = wedge_labels(1,round(size(IT_,2)/2));
        label_right = wedge_labels(end,round(size(IT_,2)/2));
        
        if label_top ~= 0
            ITcg(wedge_labels == label_top) = NaN;
        end
        if label_bot ~= 0
            ITcg(wedge_labels == label_bot) = NaN;
        end
        if label_left ~= 0
            ITcg(wedge_labels == label_left) = 0;
        end
        if label_right ~= 0
            ITcg(wedge_labels == label_right) = 0;
        end
        
        bot_top_wedge_finder = ITcg(round(size(ITcg,1)/2),:);
        
        
        ITcgc = ITcg;%IT_processor(ITcg,0,0);
        ITcgm = ITcgc - nanmean(ITcgc(:));
        IT10 = ITcgm;
        IT10(isnan(IT10)==0) = 1;
        IT10(isnan(IT10)) = 0;
        
        if wannasee_ABBAs == 1       
            AAA_i = ceil(ABBA_IT_i/length(B_range));
            BBB_i = ABBA_IT_i - (AAA_i-1)*length(B_range);

            show_IT = ITcgm;%((padding/2)+1:end-(padding/2),(padding/2)+1:end-(padding/2));
            show_IT(isnan(show_IT)) = 100;
            subplot(length(A_range),length(B_range),...
            (length(A_range)-AAA_i)*length(B_range) + BBB_i)
            imagesc(show_IT), axis xy, axis off;
            pause(.001)
        end
        
        ABBA_ITs(:,:,ABBA_IT_i) = ITcgm;
        ABBA_IT10s(:,:,ABBA_IT_i) = IT10;
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
    scone_bfc = scone_bf(b_k_1:b_k_2, b_E_1:b_E_2);
    fass = scone_bfc; 
    fass_g = mat2gray(fass);
    %fass_gm = fass_g - mean(fass_g(:));
    f_image = fass_g; 
    
    tempsize = size(ABBA_ITs(:,:,1));
    AB_corr_table = zeros([length(A_range),length(B_range)]);  
    corr_IT = zeros([(size(f_image,1)-tempsize(1)+1),(size(f_image,2)-tempsize(2)+1)]);

    for AB_IT_i = 1:(length(A_range)*length(B_range))
        AA_i = ceil(AB_IT_i / length(B_range));
        BB_i = AB_IT_i - (AA_i-1)*length(B_range);
        if ismember(AB_IT_i, nan_ABs) == 1      
            AB_corr_table(AA_i,BB_i) = NaN;
            continue
        end
        
        IT_template = ABBA_ITs(:,:,AB_IT_i);
        IT10_temp = ABBA_IT10s(:,:,AB_IT_i);
        N = size(IT_template,1);
        M = size(IT_template,2);

        corr = zeros(size(corr_IT));
        
        %{
        f_denom = zeros(size(normcorr));
        f_mean = zeros(size(normcorr));
        t_mean = nanmean(IT_template(:));
        t_denom = sum(nansum((IT_template-t_mean).^2));
        
        
        for x = 1:(size(f_image,2) - size(IT_template,2) +1)
            for y = 1:(size(f_image,1) - size(IT_template,1) +1)
                f_mean(y,x) = mean( sum( dot( IT10_temp, f_image(y:y+N-1, x:x+M-1))));
                %f_denom = sum(nansum((f_image(y:y+N-1,x:x+M-1)-f_mean).^2));
                %NC = nansum( dot( (IT_template(1:N,1:M) - t_mean),(f_image(y:y+N-1,x:x+M-1)-f_mean)));
                %normcorr(y,x) = sum(NC) / sqrt(t_denom * f_denom);
                
                
                
                
                for m = 1:size(IT_template,2)
                    %normcol = (IT_template(1:N,m) - t_mean) .* (f_image(y:y+N-1,x+m-1) - f_mean);
                    %normcorr(y,x) = normcorr(y,x) + nansum(normcol);
                    %f_denomcol = (f_image(y:y+N-1,x+m-1)-f_mean).^2;
                    %f_denom = f_denom + nansum(f_denomcol);
                    
                    for n = 1:size(IT_template,1)
                        if isnan(IT_template(n,m))==1
                            continue
                        end
                        normcorr(y,x) = normcorr(y,x) + (IT_template(n,m)-t_mean)*(f_image(y+n-1,x+m-1)-f_mean(y,x));
                        f_denom(y,x) = f_denom(y,x) + (f_image(y+n-1,x+m-1)-f_mean(y,x)).^2;
                    end
                    
                end
                %if (x == 5) && (y==5)
                %    disp(['AB_IT_i=',num2str(AB_IT_i),'t_denom=',num2str(t_denom)])
                %end
                normcorr(y,x) = normcorr(y,x) / sqrt(t_denom * f_denom(y,x));
                
            end
        end
        %}
        for x = 1:(size(f_image,2) - size(IT_template,2) +1)
            for y = 1:(size(f_image,1) - size(IT_template,1) +1)  
                f_image_under = f_image(y:y+N-1, x:x+M-1);
                corr(y,x) =  sum(nansum( IT_template .* f_image_under ));
            end
        end
        [corrmax_row,corrmax_col] = find(corr==max(corr(:)));
        AB_corr_table(AA_i,BB_i) = max(corr(:));
        
        if AB_corr_table(AA_i,BB_i) == max(AB_corr_table(:))
            corr_IT = corr;
            AA_i_it = AA_i;
            BB_i_it = BB_i;
            AB_IT_i_it = AB_IT_i;
        end 
    end
    
    
    [corrmax_row_it,corrmax_col_it] = find(corr_IT==max(corr_IT(:)));
    corr_ITg = (corr_IT - min(corr_IT(:)))/(max(corr_IT(:))-min(corr_IT(:)));
    corrspread = length(find(corr_ITg > corrspread_TH));

    DPE_b = E_0 + corrmax_col_it + b_E_1 - 1;
    DPK_b = K_0 + corrmax_row_it + b_k_1 - 1;
    DPE_pix = (DPE_b * bin_E) + bin_E/2;
    DPK_pix = (DPK_b* bin_k) + bin_k/2;

    Es_ato(i) = DPE_pix;
    ks_ato(i) = DPK_pix;
    As_ato(i) = A_range(AA_i_it);
    Bs_ato(i) = B_range(BB_i_it);
    corrs_ato(i) = max(corr_IT(:));
    corrspreads_ato(i) = corrspread;



    if wannasee == 1
        %overlay = zeros(size(f_image));
        overlay = zeros(size(scone_bf));
        overlay((1:size(IT_template,1))+corrmax_row_it-1+b_k_1-1,(1:size(IT_template,2))+corrmax_col_it-1+b_E_1-1) = ABBA_ITs(:,:,AB_IT_i_it);
        overlay(isnan(overlay)) = -.1;
        
        figure, 
        subplot(2,2,1), 
        imagesc(B_range,A_range,AB_corr_table), axis xy, hold on %title('AB\_corr\_table'), 
        plot(B_range(BB_i_it),A_range(AA_i_it),'w+'), hold off
        title(['A=',num2str(A_range(AA_i_it)),' B=',num2str(B_range(BB_i_it))])

        subplot(2,2,2),
        imagesc((1:size(corr_IT,1)),fliplr(1:size(corr_IT,2)),rot90(corr_IT,1)), axis xy, hold on
        plot(corrmax_row_it,corrmax_col_it,'w+'), hold off
        title({['spread=',num2str(corrspreads_ato(i))];['maxC=',num2str(corrs_ato(i))]})

        subplot(2,2,3), 
        imagesc((1:size(f_image,1)),fliplr(1:size(f_image,2)),rot90(f_image,1)), axis xy
        imagesc(scone_bf'), axis xy;
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




