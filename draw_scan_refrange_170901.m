tic;
scan_is = 151;%769
wannasee = 1;

bin_E = 5;
bin_k = 2;
fass_sigma = 2;

LineWidth = 1;
draw_sigma = 1;

draw_box = zeros(51,55);
E_0 = 25;%20;
K_0 = round(size(draw_box,1)/2);
draw_x = (1:size(draw_box,1))';
draw_x0 = draw_x - K_0;

A_range = 1.5:.1:2.5;%[0.5:.5:3.0];%[.6:.2:2];%[1.8];%[1,2,3];
B_range = 0:.005:.025;%[0:.002:.01];%[0,.05,.1];
E_rough_range = round(20:2:55);%round([30:2.5:60]);%round([26:2.5:50]);%round([25:2.5:50]);%[20:5:50];
K_rough_range = round(101/2 - round(size(draw_box,1)/2)) + (-20:2:20);%round([10:2:30]);%round([10:2.5:25]);%[20:5:50];
E_ref_range = -20:1:20;
K_ref_range = -10:1:10;

good_spot_A = 1.8;  %from raw_full_cone
good_spot_B = 0.025; %from raw_full_cone


A_range_eVA = A_range * (bin_E/bin_k) * .8107;
B_range_eVA = B_range * (bin_E/bin_k) * .8107;

MC_TH = 0.65;
multi_TH = 0.7;


draw_Es_904 = zeros(1,num_scans);
draw_ks_904 = zeros(1,num_scans);
draw_As_904 = zeros(1,num_scans);
draw_Bs_904 = zeros(1,num_scans);
draw_MCs_904 = zeros(1,num_scans);
draw_MCSs_904 = zeros(1,num_scans);
ABBA_rats_904 = zeros(1,num_scans);
draw_peaks_904 = zeros(1,num_scans);

ABBA_ITs = zeros(size(draw_box,1),size(draw_box,2),length(A_range)*length(B_range));






%First create all the templates%%%
ABBA_IT_i = 1;
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

        ITP = insertShape(draw_box, 'Line', draw_itp_it, 'LineWidth',LineWidth);
        if nnz(ITP(:,1))==0   %Only consider A,B values that give lines that extend across draw_box
            %ABBA_norm_table(A_i,B_i) = NaN;
            ABBA_IT_i = ABBA_IT_i + 1;
            continue
        end
        ITN = insertShape(draw_box, 'Line', draw_itn_it, 'LineWidth',LineWidth);

        ITP = mat2gray(ITP(:,:,1));
        ITN = mat2gray(ITN(:,:,1));

        IT = IT_processor(ITP + ITN, draw_sigma);   
        
        ABBA_ITs(:,:,ABBA_IT_i) = IT;
        ABBA_IT_i = ABBA_IT_i+1;
    end
end



for i = scan_is
    if rem(i,50) == 0
        disp(['Starting scan #',num2str(i)]), toc
    end
    
    cone = cones(:,:,i);
    [fass,fass_k_off] = kLOSfinder5(cone,bin_E,bin_k);
    fass = imgaussfilt(fass,fass_sigma);
    
    %%Scan first round to find where to scan more closely

    good_spot_table = zeros(size(length(K_rough_range),length(E_rough_range)));
    for E_rough_i = 1:length(E_rough_range)
        E_rough = E_rough_range(E_rough_i);

        for K_rough_i = 1:length(K_rough_range)
            K_rough = K_rough_range(K_rough_i);
            
            draww_yp = good_spot_A*(draw_x0) + good_spot_B*((draw_x0).^2) + E_0;
            draww_yn = -good_spot_A*(draw_x0) + good_spot_B*((draw_x0).^2) + E_0;

            draww_itp = horzcat(draww_yp,draw_x);
            draww_itn = horzcat(draww_yn,draw_x);

            draww_itp_curt = draww_itp(max(1,round(-good_spot_A/(2*good_spot_B)+K_0)):end,:);
            draww_itn_curt = draww_itn(1:min(length(draw_x),round(K_0+good_spot_A/(2*good_spot_B))),:);

            draww_itp_it = reshape(draww_itp_curt',1,[]);
            draww_itn_it = reshape(draww_itn_curt',1,[]);
                   
            IITP = insertShape(draw_box, 'Line', draww_itp_it, 'LineWidth',LineWidth);
            IITN = insertShape(draw_box, 'Line', draww_itn_it, 'LineWidth',LineWidth);

            IITP = mat2gray(IITP(:,:,1));
            IITN = mat2gray(IITN(:,:,1));

            IIT = IT_processor(IITP + IITN, draw_sigma);
               
            draww_scan_window = fass([1:size(draw_box,1)]+K_rough,...
                                    [1:size(draw_box,2)]+E_rough);
            
            draww_scan_window_norm = window_processor(draww_scan_window); 
           
            good_spot_table(K_rough_i,E_rough_i) = sum(dot(draww_scan_window_norm,IIT_cropped));
        end
    end
    %good_spot_table = imgaussfilt(good_spot_table,1);
    [K_scan_center_i,E_scan_center_i] = find(good_spot_table==max(good_spot_table(:)));
    %%%%%%%%%%
    K_off_range = K_rough_range(K_scan_center_i) + K_ref_range;
    E_off_range = E_rough_range(E_scan_center_i) + E_ref_range;
    %%%%%%%%%%%        
    
    
    EK_ABBAn_table = zeros(size(length(K_off_range),length(E_off_range)));
    EK_ABBAn_A_table = zeros(size(length(K_off_range),length(E_off_range)));
    EK_ABBAn_B_table = zeros(size(length(K_off_range),length(E_off_range)));

    %%%Now the actual scanning across E and k
    for E_off_i = 1:length(E_off_range)
        E_off = E_off_range(E_off_i);

        for K_off_i = 1:length(K_off_range)
            K_off = K_off_range(K_off_i);

            AB_MC_table = zeros(size(length(A_range),length(B_range)));
           
            ABBA_norm_table = zeros(size(length(A_range),length(B_range)));
            ABBA_IT_i = 1;
            for A_i = 1:length(A_range)
                A = A_range(A_i);

                for B_i = 1:length(B_range)
                    B = B_range(B_i);
                    
                    draw_scan_window = fass([1:size(draw_box,1)]+K_off,...
                                            [1:size(draw_box,2)]+E_off);
                    draw_scan_window_norm = window_processor(draw_scan_window);
                    
                    ABBA_norm_table(A_i,B_i) = sum(dot(draw_scan_window_norm,ABBA_ITs(:,:,ABBA_IT_i)));
                    ABBA_IT_i = ABBA_IT_i+1;
                end
            end 
            
            [ABBAn_it_i,BAABn_it_i] = find(ABBA_norm_table==max(ABBA_norm_table(:)),1);
            ABBAn_it = A_range(ABBAn_it_i);
            BAABn_it = B_range(BAABn_it_i);

            EK_ABBAn_table(K_off_i,E_off_i) = max(ABBA_norm_table(:));
            EK_ABBAn_A_table(K_off_i,E_off_i) = ABBAn_it;
            EK_ABBAn_B_table(K_off_i,E_off_i) = BAABn_it;      

        end    
    end

    [K_off_ABBAn_it_i, E_off_ABBAn_it_i] = find(EK_ABBAn_table==max(EK_ABBAn_table(:)));
    K_off_ABBAn_it = K_off_range(K_off_ABBAn_it_i);
    E_off_ABBAn_it = E_off_range(E_off_ABBAn_it_i);
    ABBAn_it_it = EK_ABBAn_A_table(K_off_ABBAn_it_i,E_off_ABBAn_it_i);
    BAABn_it_it = EK_ABBAn_B_table(K_off_ABBAn_it_i,E_off_ABBAn_it_i);
    MC_ABBAn_it = EK_ABBAn_table(K_off_ABBAn_it_i,E_off_ABBAn_it_i);
    
    EK_ABBAnn_table = (EK_ABBAn_table-min(EK_ABBAn_table(:)))./(max(EK_ABBAn_table(:))-min(EK_ABBAn_table(:)));
    [max_a,max_b] = find(EK_ABBAnn_table>=MC_TH);
    corr_spread = length(max_a);
    
    merp = sum(EK_ABBAnn_table(round(length(K_off_range)/2)-1:round(length(K_off_range)/2)+1,:));
    [pks,locs,wdths,proms] = findpeaks(merp,'MinPeakDistance',5);
    %figure, plot(merp)
    
    draw_Es_904(i) = (E_off_ABBAn_it + E_0) * bin_E + round(bin_E/2);
    draw_ks_904(i) = (K_off_ABBAn_it + K_0 + fass_k_off) * bin_k + round(bin_k/2);
    draw_As_904(i) = ABBAn_it_it;
    draw_Bs_904(i) = BAABn_it_it;
    draw_MCs_904(i) = MC_ABBAn_it;
    draw_MCSs_904(i) = corr_spread;
    draw_peaks_904(i) = length(pks);
    
    if wannasee == 0.5 || 1
        %%%Re-make maxcorr for A vs B scan at the chosen E_off_it,K_off_it%%%
        ABBA_norm_table_it = zeros(size(length(A_range),length(B_range)));
        ABBA_IT_i = 1;
        for A_i = 1:length(A_range)
            A = A_range(A_i);

            for B_i = 1:length(B_range)
                B = B_range(B_i);
                
                draw_scan_window_itt = fass([1:size(draw_box,1)]+K_off_ABBAn_it,...
                                [1:size(draw_box,2)]+E_off_ABBAn_it);
                draw_scan_window_norm_itt = window_processor(draw_scan_window_itt);

                ABBA_norm_table_it(A_i,B_i) = sum(dot(draw_scan_window_norm_itt,ABBA_ITs(:,:,ABBA_IT_i)));
                ABBA_IT_i = ABBA_IT_i+1;
            end
        end
        
        ABBA_split = zeros(size(ABBA_norm_table_it));
        for ABBA_i = 1:length(B_range)
            abba_i = round((ABBA_i-1)*(length(A_range)/length(B_range)))+1;
            ABBA_split(abba_i:end,ABBA_i) = 1;
        end
        
        ABBA_halfA_sum = sum(dot(ABBA_split,ABBA_norm_table_it));
        ABBA_halfB_sum = sum(dot((1-ABBA_split),ABBA_norm_table_it));
        
        ABBA_rats_904(i) = ABBA_halfA_sum / ABBA_halfB_sum;     
    end
    
    
    if wannasee == 1

        %%%%Draw Figure Showing%%%
        draw_yp_it_it = ABBAn_it_it*(draw_x-K_0) + BAABn_it_it*((draw_x-K_0).^2) + E_0;
        draw_yn_it_it = -ABBAn_it_it*(draw_x-K_0) + BAABn_it_it*((draw_x-K_0).^2) + E_0;

        draw_ittp = horzcat(draw_yp_it_it,draw_x);
        draw_ittn = horzcat(draw_yn_it_it,draw_x);


        draw_ittp_curt = draw_ittp(max(1,round(-ABBAn_it_it/(2*BAABn_it_it)+K_0)):end,:);
        draw_ittn_curt = draw_ittn(1:min(length(draw_x),round(K_0+ABBAn_it_it/(2*BAABn_it_it))),:);

        draw_ittp_it = reshape(draw_ittp_curt',1,[]);
        draw_ittn_it = reshape(draw_ittn_curt',1,[]);

        ITTP = insertShape(draw_box, 'Line', draw_ittp_it, 'LineWidth',LineWidth);
        %{
        if nnz(ITTP(:,1))==0   %Only consider A,B values that give lines that extend across draw_box
            ABBA_norm_table(A_i,B_i) = NaN;
            continue
        end
        %}
        ITTN = insertShape(draw_box, 'Line', draw_ittn_it, 'LineWidth',LineWidth);

        ITTP = mat2gray(ITTP(:,:,1));
        ITTN = mat2gray(ITTN(:,:,1));

        ITT = IT_processor(ITTP + ITTN, draw_sigma);

        draw_scan_window_it = fass([1:size(draw_box,1)]+K_off_ABBAn_it,...
                                [1:size(draw_box,2)]+E_off_ABBAn_it);
        draw_scan_window_norm_it = window_processor(draw_scan_window_it);%draw_scan_window_it ./ sum(draw_scan_window_it(:));

        draw_pic_it = zeros(size(fass));
        draw_pic_it([1:size(draw_box,1)] + K_off_ABBAn_it,...
                 [1:size(draw_box,2)] + E_off_ABBAn_it ) = ITT;
        
             
        K_conf_range = K_off_ABBAn_it + [-10:1:10];
        E_conf_range = E_off_ABBAn_it + [-20:1:20];
        find_multi_ks_table = zeros(length(K_conf_range),length(E_conf_range));
        for E_conf_i = 1:length(E_conf_range)
            E_conf = E_conf_range(E_conf_i);
            for K_conf_i = 1:length(K_conf_range)
                K_conf = K_conf_range(K_conf_i);
                
                drawww_scan_window = fass([1:size(draw_box,1)]+K_conf,...
                                    [1:size(draw_box,2)]+E_conf);
                drawww_scan_window_norm = window_processor(drawww_scan_window);%drawww_scan_window ./ sum(drawww_scan_window(:));
              
                ABBA_IT_i = find(A_range == ABBAn_it_it);
                BAAB_IT_i = find(B_range == BAABn_it_it);
                find_multi_ks_table(K_conf_i,E_conf_i) = sum(dot(drawww_scan_window_norm,ABBA_ITs(:,:,(ABBA_IT_i-1)*length(B_range)+BAAB_IT_i)));
            end
        end
        
        gray_multi_ks_table = mat2gray(find_multi_ks_table);
        ASDF = length(find(gray_multi_ks_table >= multi_TH));
        
        figure,
        
        ax5 = subplot(3,2,5);
        imagesc(draw_scan_window_it), axis xy
        colormap(ax5, jet)
        title('Best Window')
        
        %ax6 = subplot(3,2,6);
        %imagesc(ITT_cropped), axis xy
        %colormap(ax6, bone)
        %title('Best Draw Fit')
        
        ax2 = subplot(3,2,2);
        imagesc(gray_multi_ks_table), axis xy
        colormap(ax2, jet)
        title(['Find multi, above ',num2str(multi_TH),'=',num2str(ASDF)])
        
        ax1 = subplot(3,2,1);
        FixedWidth = get(0,'FixedWidthFontName');
        %text(0.2,0.5,{['A (eVA) :  ',num2str(A_it_it*bin_E/bin_k*.8107)];['B (eVA) :  ',num2str(B_it_it*bin_E/bin_k*.8107)];['MC:  ',num2str(MC_it)];...
        text(0,0.5,{...
            ['A (bpix units):          ',num2str(ABBAn_it_it)];...
            ['B (bpix units):          ',num2str(BAABn_it_it)];...
            ['MC (faswin normed):      ',num2str(MC_ABBAn_it)];...
            ['     '];...
            ['DP E coor (bpix units):  ',num2str(E_off_ABBAn_it+E_0)];...
            ['DP K coor (bpix units):  ',num2str(K_off_ABBAn_it+K_0)]} ,'FontName',FixedWidth,'FontSize',9);
        axis off

        ax3 = subplot(3,2,3); 
        imagesc(EK_ABBAn_table), axis xy
        colormap(ax3, jet)
        title(['Max Corr, above ',num2str(MC_TH),' = ',num2str(corr_spread)])
        %yticks([1,length(K_off_range)])
        %yticklabels({num2str(K_off_range(1)+K_0),num2str(K_off_range(end)+K_0)})
        %ylabel('DP K (binned pixels, fass frame)','FontSize',8)
        %xticks([1,length(E_off_range)])
        %xticklabels({num2str(E_off_range(1)+E_0),num2str(E_off_range(end)+E_0)});
        %xlabel('DP E (binned pixels, fass frame)','FontSize',8)

        ax4 = subplot(3,2,4);
        imagesc(ABBA_norm_table_it), axis xy, title('MaxCorr for A-vs-B at DP')
        colormap(ax4, jet)
       
        xticks([1:2:length(B_range)])
        xticklabels(num2str(B_range(1:2:length(B_range))))
        
        set(gca, 'XTick',[1:2:length(B_range)], 'XTickLabel',{num2str(B_range_eVA(1,1:2:length(B_range))')}, 'FontSize',6)
        set(gca, 'YTick',[1:2:length(A_range)], 'YTickLabel',{num2str(A_range_eVA(1,1:2:length(A_range))')}, 'FontSize',6)
        set(gca, 'XTick',[1:2:length(B_range)], 'XTickLabel',{num2str(B_range(1,1:2:length(B_range))')}, 'FontSize',6)
        set(gca, 'YTick',[1:2:length(A_range)], 'YTickLabel',{num2str(A_range(1,1:2:length(A_range))')}, 'FontSize',6)
        yticks([1:2:length(A_range)])
        A_range_ = A_range';
        yticklabels(num2str(A_range_(1:2:length(A_range))))
        xlabel('B term (pixel units)','FontSize',8)
        ylabel('A term (pixel units)','FontSize',8)

        ax6 = subplot(3,2,6);
        E = mat2gray(fass(16:end-15,:));
        I = mat2gray(draw_pic_it(16:end-15,:));
        greenlol = cat(3,ones(size(E)),zeros(size(E)),zeros(size(E)));
        imshow(E,'InitialMag','fit','DisplayRange',[0 .5]), hold on;
        h = imshow(greenlol); axis xy 
        set(h,'AlphaData',I);
        hold off;
        colormap(ax6, winter)
        
        suptitle(['Scan i=',num2str(i)])
        pause(.01)
    end
end

toc;

function processed_IT = IT_processor(input_IT, draw_sigma)
processed_IT_ = imgaussfilt(input_IT,draw_sigma);
processed_IT_ = mat2gray(processed_IT_);
processed_IT = processed_IT_ ./ sum(processed_IT_(:));
end

function processed_window = window_processor(input_window)
processed_window = input_window ./ sum(input_window(:));
end


