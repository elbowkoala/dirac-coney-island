tic;
%%%%%PARAMS%%%%%%%%%%%%%%%%
scan_is =  1:num_scans;%round(961*rand);% [round(961*rand), round(961*rand), round(961*rand), round(961*rand)];
wannasee = 0;

bin_E = 5;
bin_k = 2;
fass_sigma = 3;
b_E_1 = 25;
b_E_2 = 130;
b_k_1 = 40;
b_k_2 = 140;

IT_sigma = 3;

LineWidth = 1;
draw_sigma = 4;

draw_box = zeros(71,75);
E_0 = 30;
K_0 = round(size(draw_box,1)/2);
draw_x = (1:size(draw_box,1))';
draw_x0 = draw_x - K_0;

A_range = [1.5:.1:2.5];%[0.5:.5:3.0];%[.6:.2:2];%[1.8];%[1,2,3];
B_range = [0:.002:.02];%[0:.002:.01];%[0,.05,.1];

A_range_eVA = A_range * (bin_E/bin_k) * .8107;
B_range_eVA = B_range * (bin_E/bin_k) * .8107;

corrspread_TH = 0.1;
%%%%%%%%%%%%%%%%%%%%

dssd_Es = zeros(1,num_scans);
dssd_ks = zeros(1,num_scans);
dssd_As = zeros(1,num_scans);
dssd_Bs = zeros(1,num_scans);
dssd_corrs = zeros(1,num_scans);
dssd_corrspreads = zeros(1,num_scans);


%First draw all the A,B templates%%%
ABBA_ITs = zeros(size(draw_box,1),size(draw_box,2),length(A_range)*length(B_range));
nan_ABs = [];
ABBA_IT_i = 1;
nan_AB_i = 1;
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
        
        IT = imgaussfilt( ITP+ITN, IT_sigma);
        ITg = mat2gray(IT);
        ITgm = ITg - mean(ITg(:));
        
        IT_ = IT;
        IT_((IT~=0))=1;
        IT_ = abs(1-IT_);
        
        IT_wedge = bwconncomp(IT_);
        if IT_wedge.NumObjects == 3
            ITgm(IT_wedge.PixelIdxList{3}) = NaN;
        elseif IT_wedge.NumObjects == 4
            ITgm(IT_wedge.PixelIdxList{2}) = NaN;
            ITgm(IT_wedge.PixelIdxList{4}) = NaN;
        end
        
        AAA_i = ceil(ABBA_IT_i/length(B_range));
        BBB_i = ABBA_IT_i - (AAA_i-1)*length(B_range);
       
        %subplot(length(A_range),length(B_range),...
        %(length(A_range)-AAA_i)*length(B_range) + BBB_i)
        %imshow(IT_template,'InitialMagnification','fit'), axis xy
        %pause(.001)

        ABBA_ITs(:,:,ABBA_IT_i) = ITgm;
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
    fass_gm = fass_g - mean(fass_g(:));
    fass_image = fass_gm; 
    
    AB_corr_table = zeros([length(A_range),length(B_range)]);  
    for AB_IT_i = 1:(length(A_range)*length(B_range))
        AA_i = ceil(AB_IT_i / length(B_range));
        BB_i = AB_IT_i - (AA_i-1)*length(B_range);
        if ismember(AB_IT_i, nan_ABs) == 1      
            AB_corr_table(AA_i,BB_i) = NaN;
            continue
        end
        
        IT_template = ABBA_ITs(:,:,AB_IT_i);

        corr = zeros([(size(fass_image,1)-size(IT_template,1)+1),(size(fass_image,2)-size(IT_template,2)+1)]);

        for x = 1:(size(fass_image,2) - size(IT_template,2) +1)
            for y = 1:(size(fass_image,1) - size(IT_template,1) +1)
                for m = 1:size(IT_template,2)
                    for n = 1:size(IT_template,1)
                        if isnan(IT_template(n,m))==1
                            continue
                        end
                        corr(y,x) = corr(y,x) + (fass_image(y+n-1, x+m-1) - IT_template(n,m)).^2;
                    end
                end
            end
        end

        [corrmin_row,corrmin_col] = find(corr==min(corr(:)));
        AB_corr_table(AA_i,BB_i) = min(corr(:));
        
        if AB_corr_table(AA_i,BB_i) == min(AB_corr_table(AB_corr_table~=0))
            corr_IT = corr;
            AA_it_i = AA_i;
            BB_it_i = BB_i;
            AB_IT_i_it = AB_IT_i;
        end
        
    end
    
    [ITmin_row,ITmin_col] = find(corr_IT==min(corr_IT(:)));
    corr_ITg = (corr_IT - min(corr_IT(:)))/(max(corr_IT(:))-min(corr_IT(:)));
    corrspread = length(find(corr_ITg < corrspread_TH));
   
    DPE_b = E_0 + ITmin_col + b_E_1 - 1;
    DPK_b = K_0 + ITmin_row + b_k_1 - 1;
    DPE_pix = (DPE_b * bin_E) + bin_E/2;
    DPK_pix = (DPK_b* bin_k) + bin_k/2;
    
    dssd_Es(i) = DPE_pix;
    dssd_ks(i) = DPK_pix;
    dssd_As(i) = A_range(AA_it_i);
    dssd_Bs(i) = B_range(BB_it_i);
    dssd_corrs(i) = min(corr_IT(:));
    dssd_corrspreads(i) = corrspread;
    

    
    if wannasee == 1
        overlay = zeros(size(fass_image));
        overlay((1:size(IT_template,1))+ITmin_row-1,(1:size(IT_template,2))+ITmin_col-1) = ABBA_ITs(:,:,AB_IT_i_it);
        overlay(isnan(overlay)) = 0;
        
        figure, 
        subplot(2,2,1), 
        imagesc(B_range,A_range,AB_corr_table), axis xy, title('AB\_corr\_table'), hold on;
        plot(B_range(BB_it_i),A_range(AA_it_i),'w+'), hold off

        subplot(2,2,2),
        imagesc((1:size(corr_IT,1)),fliplr(1:size(corr_IT,2)),rot90(corr_IT,1)), axis xy, hold on
        plot(ITmin_row,ITmin_col,'w+'), hold off
        title({['spread=',num2str(dssd_corrspreads(i))];['minC=',num2str(min(corr_IT(:)))]})

        subplot(2,2,3), 
        imagesc((1:size(fass_image,1)),fliplr(1:size(fass_image,2)),rot90(fass_image,1)), axis xy
        title(['DPI=',num2str(DPI_big(i))])

        subplot(2,2,4), imagesc((1:size(fass_image,1)),fliplr(1:size(fass_image,2)),rot90(fass_image+1.5*overlay,1)), axis xy, hold on;
        plot([1,size(IT_template,1),size(IT_template,1),1,1]+ITmin_row-1,[1,1,size(IT_template,2),size(IT_template,2),1]+ITmin_col-1,'w'), hold off;
        title(['Epix=',num2str(DPE_pix),' Kpix=',num2str(DPK_pix)])
        suptitle(['i=',num2str(i)])
    end
    
end

toc;




