tic;
scan_is =  1;%round(rand*961);
%769
wannasee = 1;

bin_E = 5;
bin_k = 2;
fass_sigma = 2.5;

LineWidth = 1;
draw_sigma = 2.5;
krillin_sigma = 1;

draw_box = zeros(61,61);
E_0 = 20;%20;
E_B1 =2;
E_B2 = size(draw_box,2)-3;
K_0 = round(size(draw_box,1)/2);
draw_x = (1:size(draw_box,1))';
draw_x0 = draw_x - K_0;

A_range = [1.5:.05:2.1];%[0.5:.5:3.0];%[.6:.2:2];%[1.8];%[1,2,3];
B_range = 0:.001:.005;%[0:.002:.01];%[0,.05,.1];
E_rough_range = round(35:2:60);
K_rough_range = round(101/2 - round(size(draw_box,1)/2)) + (-20:2:20);
E_ref_range = -8:1:8;
K_ref_range = -3:1:3;
E_conff_range = -5:1:5;
K_conff_range = -5:1:5;

rough_scan_A = 1.85; 
rough_scan_B = 0.0; 
%from raw_full_cone 9/5/17, use rough_scan_A=1.9,B=0

A_range_eVA = A_range * (bin_E/bin_k) * .8107;
B_range_eVA = B_range * (bin_E/bin_k) * .8107;

MC_TH = 0.65;
multi_TH = 0.7;

draw_Es_915d = zeros(1,num_scans);
draw_ks_915d = zeros(1,num_scans);
draw_As_915d = zeros(1,num_scans);
draw_Bs_915d = zeros(1,num_scans);
draw_MCs_915d = zeros(1,num_scans);
draw_MCSs_915d = zeros(1,num_scans);
draw_rats_915d = zeros(1,num_scans);
piccolos_915d = zeros(1,num_scans);
%}
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
            %ABBA_norm_table(A_i,B_i) = NaN;
            ABBA_IT_i = ABBA_IT_i + 1;
            continue
        end
        ITN = insertShape(draw_box, 'Line', draw_itn_it, 'LineWidth',LineWidth);

        ITP = mat2gray(ITP(:,:,1));
        ITN = mat2gray(ITN(:,:,1));

        IT = IT_processor(ITP + ITN, draw_sigma, E_B1,E_0,E_B2);   
        
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
    rough_scan_table = zeros(size(length(K_rough_range),length(E_rough_range)));
    rough_scan_table_ratios = zeros(size(length(K_rough_range),length(E_rough_range)));  
    for E_rough_i = 1:length(E_rough_range)
        E_rough = E_rough_range(E_rough_i);

        for K_rough_i = 1:length(K_rough_range)
            K_rough = K_rough_range(K_rough_i);
            
            rough_scan_A_i = find(A_range == rough_scan_A);
            rough_scan_B_i = find(B_range == rough_scan_B);
            rough_scan_IT = ABBA_ITs(:,:,(rough_scan_A_i-1)*length(B_range)+rough_scan_B_i);
            
            draww_scan_window = fass([1:size(draw_box,1)]+K_rough,...
                                    [1:size(draw_box,2)]+E_rough);
            
            draww_scan_window_norm = window_processor(draww_scan_window); 
           
            rough_scan_table(K_rough_i,E_rough_i) = sum(dot(draww_scan_window_norm,rough_scan_IT));
            
            fernbot = sum(draww_scan_window_norm(1:K_0,:));
            ferntop = sum(draww_scan_window_norm(K_0+1:end,:));          
            
            rough_scan_table(K_rough_i,E_rough_i) = sum(dot(draww_scan_window_norm,rough_scan_IT));
            rough_scan_table_ratios(K_rough_i,E_rough_i) = mean(ferntop./fernbot);

        end
    end
    %{
    figure, 
    subplot(2,4,1), imagesc(rough_scan_table), axis xy
    title('total dot products')
    xchikkusu = 1:3:length(E_rough_range);
    ychikkusu = 1:2:length(K_rough_range);
    xraberus = cell(1,length(xchikkusu));
    yraberus = cell(1,length(ychikkusu));
    for i = 1:length(xchikkusu)
        xraberus{i} = num2str(E_rough_range(xchikkusu(i)));
    end
    for i = 1:length(ychikkusu)
        yraberus{i} = num2str(K_rough_range(ychikkusu(i)));
    end
    xticks(xchikkusu)
    yticks(ychikkusu)
    xticklabels(xraberus)
    yticklabels(yraberus)
    xlabel('E pix')
    ylabel('K pix')
    
    subplot(2,4,2), imagesc(abs(1-rough_scan_table_ratios)), axis xy
    title('top/bot ratios')
    xchikkusu = 1:3:length(E_rough_range);
    ychikkusu = 1:2:length(K_rough_range);
    xraberus = cell(1,length(xchikkusu));
    yraberus = cell(1,length(ychikkusu));
    for i = 1:length(xchikkusu)
        xraberus{i} = num2str(E_rough_range(xchikkusu(i)));
    end
    for i = 1:length(ychikkusu)
        yraberus{i} = num2str(K_rough_range(ychikkusu(i)));
    end
    xticks(xchikkusu)
    yticks(ychikkusu)
    xticklabels(xraberus)
    yticklabels(yraberus)
    xlabel('E pix')
    ylabel('K pix')
    
    subplot(2,4,3), imagesc(rough_scan_table_bothalf), axis xy
    title('bot half dot prods')
    xchikkusu = 1:3:length(E_rough_range);
    ychikkusu = 1:2:length(K_rough_range);
    xraberus = cell(1,length(xchikkusu));
    yraberus = cell(1,length(ychikkusu));
    for i = 1:length(xchikkusu)
        xraberus{i} = num2str(E_rough_range(xchikkusu(i)));
    end
    for i = 1:length(ychikkusu)
        yraberus{i} = num2str(K_rough_range(ychikkusu(i)));
    end
    xticks(xchikkusu)
    yticks(ychikkusu)
    xticklabels(xraberus)
    yticklabels(yraberus)
    xlabel('E pix')
    ylabel('K pix')
    
    subplot(2,4,4), imagesc(abs(1-rough_scan_table_tophalf./rough_scan_table_bothalf)), axis xy
    title('ratio top/bot dot prods')
    xchikkusu = 1:3:length(E_rough_range);
    ychikkusu = 1:2:length(K_rough_range);
    xraberus = cell(1,length(xchikkusu));
    yraberus = cell(1,length(ychikkusu));
    for i = 1:length(xchikkusu)
        xraberus{i} = num2str(E_rough_range(xchikkusu(i)));
    end
    for i = 1:length(ychikkusu)
        yraberus{i} = num2str(K_rough_range(ychikkusu(i)));
    end
    xticks(xchikkusu)
    yticks(ychikkusu)
    xticklabels(xraberus)
    yticklabels(yraberus)
    xlabel('E pix')
    ylabel('K pix')
    
    subplot(2,3,[4 5 6])
    imagesc(fass), axis xy, hold on;
    plot(E_0+[E_rough_range(1),E_rough_range(end),E_rough_range(end),E_rough_range(1),E_rough_range(1)],...
        K_0+[K_rough_range(1),K_rough_range(1),K_rough_range(end),K_rough_range(end),K_rough_range(1)],'r'), hold off
    %}
    
    [K_scan_center_i,E_scan_center_i] = find(rough_scan_table==max(rough_scan_table(:)));
    
    K_off_range = K_rough_range(K_scan_center_i) + K_ref_range;
    E_off_range = E_rough_range(E_scan_center_i) + E_ref_range;   
    
    EK_AB_MC_table = zeros(size(length(K_off_range),length(E_off_range)));
    EK_AB_A_table = zeros(size(length(K_off_range),length(E_off_range)));
    EK_AB_B_table = zeros(size(length(K_off_range),length(E_off_range)));

    %%%Now the actual scanning across E and k
    for E_off_i = 1:length(E_off_range)
        E_off = E_off_range(E_off_i);

        for K_off_i = 1:length(K_off_range)
            K_off = K_off_range(K_off_i);
            
            %Continue to next E,K coordinate if current one makes fit box
            %out of bounds of fass
            if (1+K_off<=0) || (size(draw_box,1)+K_off>=size(fass,1)) || (1+E_off<=0) || (size(draw_box,2)+E_off>=size(fass,2));
                EK_AB_MC_table(K_off_i,E_off_i) = 0;
                EK_AB_A_table(K_off_i,E_off_i) = 0;
                EK_AB_B_table(K_off_i,E_off_i) = 0;      

            continue
                    
            %Create corr table for AvsB at each E,K point
            AB_MC_table = zeros(size(length(A_range),length(B_range)));           
            AB_MC_norm_table = zeros(size(length(A_range),length(B_range)));
            ABBA_IT_i = 1;
            for A_i = 1:length(A_range)
                A = A_range(A_i);
                for B_i = 1:length(B_range)
                    B = B_range(B_i);          
                        
                    draw_scan_window = fass([1:size(draw_box,1)]+K_off,...
                                            [1:size(draw_box,2)]+E_off);
                    draw_scan_window_norm = window_processor(draw_scan_window);
                    
                    AB_MC_norm_table(A_i,B_i) = sum(dot(draw_scan_window_norm, ABBA_ITs(:,:,ABBA_IT_i)));
                    ABBA_IT_i = ABBA_IT_i+1;
                end
            end 
            %Choose max corr A,B at the E,K point, store in table
            [A_it_i,B_it_i] = find(AB_MC_norm_table==max(AB_MC_norm_table(:)),1);
            A_it = A_range(A_it_i);
            B_it = B_range(B_it_i);

            EK_ABBAn_table(K_off_i,E_off_i) = max(AB_MC_norm_table(:));
            EK_ABBAn_A_table(K_off_i,E_off_i) = A_it;
            EK_ABBAn_B_table(K_off_i,E_off_i) = B_it;      
        end    
    end
    
    %Find E,k with max corr value, then ID which A,B values were used there
    [K_off_it_i, E_off_it_i] = find(EK_AB_MC_table==max(EK_AB_MC_table(:)));
    K_off_it = K_off_range(K_off_it_i);
    E_off_it = E_off_range(E_off_it_i);       
   
    A_it_it = EK_AB_A_table(K_off_it_i,E_off_it_i);
    B_it_it = EK_AB_B_table(K_off_it_i,E_off_it_i);
    MC_it = EK_AB_MC_table(K_off_it_i,E_off_it_i);
    
    EK_ABBAnn_table = (EK_AB_MC_table-min(EK_AB_MC_table(:)))./(max(EK_AB_MC_table(:))-min(EK_AB_MC_table(:)));
    [max_a,max_b] = find(EK_ABBAnn_table>=MC_TH);
    corr_spread = length(max_a);
   
    %Draw the chosen A,B curves and the chosen fass window
    ABBA_IT_i = find(A_range == A_it_it);
    BAAB_IT_i = find(B_range == B_it_it);
    IT_IT = ABBA_ITs(:,:,(ABBA_IT_i-1)*length(B_range)+BAAB_IT_i);
    
    fass_window_it = fass([1:size(draw_box,1)]+K_off_it,...
                            [1:size(draw_box,2)]+E_off_it);
    fass_window_norm_it = window_processor(fass_window_it); 
    
    
    krillin_yp = A_it_it*(draw_x0) + B_it_it*((draw_x0).^2) + E_0;
    krillin_yn = -A_it_it*(draw_x0) + B_it_it*((draw_x0).^2) + E_0;

    krillin_itp = horzcat(krillin_yp,draw_x0+K_0);
    krillin_itn = horzcat(krillin_yn,draw_x0+K_0);

    krillin_itp_curt = krillin_itp(max(1,round(-A_it_it/(2*B_it_it))+K_0):end,:);
    krillin_itn_curt = krillin_itn(1:min(length(draw_x0),round(K_0+A_it_it/(2*B_it_it))),:);

    krillin_itp_it = reshape(krillin_itp_curt',1,[]);
    krillin_itn_it = reshape(krillin_itn_curt',1,[]);

    krillin_ITP = insertShape(draw_box, 'Line', krillin_itp_it);
    krillin_ITN = insertShape(draw_box, 'Line', krillin_itn_it);
    krillin_ITP = krillin_ITP(:,:,1);
    krillin_ITN = krillin_ITN(:,:,1);
    
    krillin = krillin_ITP + krillin_ITN; 
    krillin = imgaussfilt(krillin, krillin_sigma); %to not penalize intensity from a good surface band that's slightly broad
    krillin(krillin~=0)=1;
    krillin = abs(1-krillin);

    gohan = bwconncomp(krillin);
    krillin(gohan.PixelIdxList{1})=0;
    krillin(gohan.PixelIdxList{2})=0;
    krillin(gohan.PixelIdxList{3})=0;
    
    piccolo = sum(dot(krillin,fass_window_it));
    
    sumbot = sum(fass_window_norm_it(1:K_0,1:E_0+10));
    sumtop = sum(fass_window_norm_it(K_0+1:end,1:E_0+10));
    
    sumtoptrash = find(sumtop==0);
    trashin = 0;
    for trashi = sumtoptrash
        sumbot(trashi-trashin) = [];
        sumtop(trashi-trashin) = [];
        trashin = trashin+1;
    end
    mean_rats = nanmean(sumbot./sumtop);
    
    %Make the AvsB scan at the chosen E,k point
    ABatDP = zeros(length(A_range), length(B_range));
    ABatDP_i = 1;
    for aa_i = 1:length(A_range)
        aa = A_range(aa_i);
        for bb = 1:length(B_range)
            bb = B_range(bb_i);
            ABatDP(aa_i,bb_i) = sum(dot(draw_window_norm_it,ABBA_ITs(:,:,ABatDP_i)));
            ABatDP_i = ABatDP_i+1;
        end
    end
    
            
    
    draw_Es_915d(i) = (E_off_it + E_0) * bin_E + round(bin_E/2);
    draw_ks_915d(i) = (K_off_it + K_0 + fass_k_off) * bin_k + round(bin_k/2);
    draw_As_915d(i) = A_it_it;
    draw_Bs_915d(i) = B_it_it;
    draw_MCs_915d(i) = MC_it;
    draw_MCSs_915d(i) = corr_spread;
    
    draw_rats_915d(i) = abs(1-mean_rats);
    piccolos_915d(i) = piccolo;
      
    if wannasee == 1
        draw_pic_it = zeros(size(fass));
        draw_pic_it([1:size(draw_box,1)] + K_off_it,...
                 [1:size(draw_box,2)] + E_off_it ) = IT_IT;
             
        K_conf_range = K_off_it + K_conff_range;
        E_conf_range = E_off_it + E_conff_range;
        smear_table = zeros(length(K_conf_range),length(E_conf_range));
        for E_conf_i = 1:length(E_conf_range)
            E_conf = E_conf_range(E_conf_i);
            for K_conf_i = 1:length(K_conf_range)
                K_conf = K_conf_range(K_conf_i);
                
                if (1+K_conf<=0) || (size(draw_box,1)+K_conf>=size(fass,1)) || (1+E_conf<=0) || (size(draw_box,2)+E_conf>=size(fass,2));
                    smear_table(K_conf_i,E_conf_i) = 0;
                    continue
                end
                drawww_scan_window = fass([1:size(draw_box,1)]+K_conf,...
                                    [1:size(draw_box,2)]+E_conf);
                drawww_scan_window_norm = window_processor(drawww_scan_window);
                
                %ABBA_IT_i = find(A_range == ABBAn_it_it);
                %BAAB_IT_i = find(B_range == BAABn_it_it);
                smear_table(K_conf_i,E_conf_i) = sum(dot(drawww_scan_window_norm,ABBA_ITs(:,:,(ABBA_IT_i-1)*length(B_range)+BAAB_IT_i)));
            end
        end
        
        gray_smear_table = mat2gray(smear_table);
        the_smear = length(find(gray_smear_table >= multi_TH));
        
        
        
        %%%%%%%Draw Figure%%%%%%%%%%%
        figure,
        
        ax5 = subplot(2,4,2);
        draw_scan_window_norm_ita = rot90(fass_window_norm_it,-1);
        imagesc(draw_scan_window_norm_ita), axis xy
        colormap(ax5, jet)
        title('Best Window','FontSize',8)
        
        %ax6 = subplot(3,2,6);
        %imagesc(ITT_cropped), axis xy
        %colormap(ax6, bone)
        %title('Best Draw Fit')
        
        ax2 = subplot(2,4,5);
        rough_scan_tablea = rot90(rough_scan_table,-1);
        %imagesc(rough_scan_tablea), axis xy
        imagesc(AB_MC_norm_table), axis xy
        xlabel('square (B) term')
        ylabel('linear (A) term')
        colormap(ax2, jet)
        %title('Rough scan','FontSize',8)
        title('A vs B at DP','FontSize',8)
        %title(['Find multi, above ',num2str(multi_TH),'=',num2str(ASDF)],'FontSize',8)
        
        ax1 = subplot(2,4,1);
        FixedWidth = get(0,'FixedWidthFontName');
        %text(0.2,0.5,{['A (eVA) :  ',num2str(A_it_it*bin_E/bin_k*.8107)];['B (eVA) :  ',num2str(B_it_it*bin_E/bin_k*.8107)];['MC:  ',num2str(MC_it)];...
        text(0,0.5,{...
            ['A(bpix):  ',num2str(A_it_it)];...
            ['B(bpix):  ',num2str(B_it_it)];...
            ['MC(arb):  ',num2str(MC_it)];...
            ['     '];...
            ['E(bpix):  ',num2str(E_off_it+E_0)];...
            ['K(bpix):  ',num2str(K_off_it+K_0)]} ,'FontName',FixedWidth,'FontSize',8);
        axis off

        ax3 = subplot(2,4,6); 
        gray_smear_tablea = rot90(gray_smear_table,-1);
        imagesc(gray_smear_tablea), axis xy
        colormap(ax3, jet)
        title({['Ref Scan'];['Pts > ',num2str(multi_TH),': ',num2str(the_smear)]},'FontSize',8)
        %yticks([1,length(K_off_range)])
        %yticklabels({num2str(K_off_range(1)+K_0),num2str(K_off_range(end)+K_0)})
        %ylabel('DP K (binned pixels, fass frame)','FontSize',8)
        %xticks([1,length(E_off_range)])
        %xticklabels({num2str(E_off_range(1)+E_0),num2str(E_off_range(end)+E_0)});
        %xlabel('DP E (binned pixels, fass frame)','FontSize',8)

        %ax4 = subplot(3,2,4);
        %imagesc(ABBA_norm_table_it), axis xy, title('MaxCorr for A-vs-B at DP')
        %colormap(ax4, jet)
       
        %xticks([1:2:length(B_range)])
        %xticklabels(num2str(B_range(1:2:length(B_range))))
        
        %set(gca, 'XTick',[1:2:length(B_range)], 'XTickLabel',{num2str(B_range_eVA(1,1:2:length(B_range))')}, 'FontSize',6)
        %set(gca, 'YTick',[1:2:length(A_range)], 'YTickLabel',{num2str(A_range_eVA(1,1:2:length(A_range))')}, 'FontSize',6)
        %set(gca, 'XTick',[1:2:length(B_range)], 'XTickLabel',{num2str(B_range(1,1:2:length(B_range))')}, 'FontSize',6)
        %set(gca, 'YTick',[1:2:length(A_range)], 'YTickLabel',{num2str(A_range(1,1:2:length(A_range))')}, 'FontSize',6)
        %yticks([1:2:length(A_range)])
        %A_range_ = A_range';
        %yticklabels(num2str(A_range_(1:2:length(A_range))))
        %xlabel('B term (pixel units)','FontSize',8)
        %ylabel('A term (pixel units)','FontSize',8)

        ax6 = subplot(2,4,[3,4,7,8]);
        E = mat2gray(norman(fass(16:end-15,:),0,30));
        Ea = rot90(E,-1);
        I = mat2gray(draw_pic_it(16:end-15,:));
        Ia = rot90(I,-1);
        greenlol = cat(3,ones(size(Ea)),zeros(size(Ea)),zeros(size(Ea)));
        imshow(Ea,'InitialMag','fit','DisplayRange',[0 .5]), hold on;
        h = imshow(greenlol);  axis xy
        set(h,'AlphaData',Ia);
        hold off;
        %colormap(ax6, winter)
        
        suptitle({['Scan i=',num2str(i)]})
        pause(.01)
    end
end

toc;




