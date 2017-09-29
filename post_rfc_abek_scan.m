tic;
%%%%%PARAMS%%%%%%%%%%%%%%%%
scan_is = [round(961*rand), round(961*rand), round(961*rand), round(961*rand)];
wannasee = 1;

bin_E = 5;
bin_k = 2;
fass_sigma = 2;

LineWidth = 1;
draw_sigma = 2;
krillin_sigma = 0.5;

draw_box = zeros(61,61);
E_0 = 20;
E_B1 =5;
E_B2 = size(draw_box,2)-10;
K_0 = round(size(draw_box,1)/2);
draw_x = (1:size(draw_box,1))';
draw_x0 = draw_x - K_0;

A_range = [1.5:.05:2.5];
B_range = 0:.001:.01;
E_rough_range = (40:100) - E_0;
K_rough_range = (80:100) - K_0;
E_ref_range = -5:1:5;
K_ref_range = -5:1:5;
E_conff_range = -10:1:10;
K_conff_range = -10:1:10;

rough_scan_A = 1.9; 
rough_scan_B = 0.0; 
if (ismember(rough_scan_A,A_range) == 0) || (ismember(rough_scan_B,B_range) == 0)
    disp('Please pick a rough starting point (A,B) that is included in [A,B_range]')
    return
end

A_range_eVA = A_range * (bin_E/bin_k) * (pix2eV/pix2invA);
B_range_eVA = B_range * (bin_E/bin_k) * (pix2eV/pix2invA);

MC_TH = 0.65;
multi_TH = 0.65;
%%%%%%%%%%%%%%%%%%%%

%{
%Create empty tables%
ABEK_Es = zeros(1,num_scans);
ABEK_ks = zeros(1,num_scans);
ABEK_As = zeros(1,num_scans);
ABEK_Bs = zeros(1,num_scans);
ABEK_MCs = zeros(1,num_scans);
ABEK_MCSs = zeros(1,num_scans);
combadges = zeros(1,num_scans);
TBDs = zeros(1,num_scans);
ABjudge = zeros(1,num_scans);
%}
ABBA_ITs = zeros(size(draw_box,1),size(draw_box,2),length(A_range)*length(B_range));

%First draw all the A,B templates%%%
A_IT_i = 1;
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
            A_IT_i = A_IT_i + 1;
            continue
        end
        ITN = insertShape(draw_box, 'Line', draw_itn_it, 'LineWidth',LineWidth);
        
        ITP = ITP(:,:,1);
        ITN = ITN(:,:,1);
        
        IT = IT_processor(ITP + ITN, draw_sigma, E_B1,E_0,E_B2);   
        
        ABBA_ITs(:,:,A_IT_i) = IT;
        A_IT_i = A_IT_i+1;
    end
end


%%%For each individual scan image.................................
for i = scan_is
    if rem(i,50) == 0
        disp(['Starting scan #',num2str(i)]), toc
    end
    
    cone = cones(:,:,i);
    coneb = Binning_2d(cone,bin_E,bin_k);
    conebf = imgaussfilt(coneb,fass_sigma); %[fass,fass_k_off] = kLOSfinder5(cone,bin_E,bin_k);
    %fass = imgaussfilt(fass,fass_sigma);
    fass = conebf;
    
    %%Scan first round roughly to find where to scan more closely, tbd is diff in
    %%summed intensities in top half of fass window vs bottom
    rough_scan_table = zeros(size(length(K_rough_range),length(E_rough_range)));
    rough_scan_table_tbd = zeros(size(length(K_rough_range),length(E_rough_range)));  
    for E_rough_i = 1:length(E_rough_range)
        E_rough = E_rough_range(E_rough_i);

        for K_rough_i = 1:length(K_rough_range)
            K_rough = K_rough_range(K_rough_i);
            
            rough_scan_A_i = find(A_range == rough_scan_A);
            rough_scan_B_i = find(B_range == rough_scan_B);
            rough_scan_IT = ABBA_ITs(:,:,(rough_scan_A_i-1)*length(B_range)+rough_scan_B_i);
            
            rough_fass_window = fass([1:size(draw_box,1)]+K_rough, [1:size(draw_box,2)]+E_rough);
            rough_fass_window_norm = window_processor(rough_fass_window); 
           
            rough_scan_table(K_rough_i,E_rough_i) = sum(dot(rough_fass_window_norm,rough_scan_IT));
            
            sumIbot = sum(rough_fass_window_norm(1:K_0,:));
            sumItop = sum(rough_fass_window_norm(K_0+1:end,:));          
            
            rough_scan_table(K_rough_i,E_rough_i) = sum(dot(rough_fass_window_norm,rough_scan_IT));
            rough_scan_table_tbd(K_rough_i,E_rough_i) = mean(abs(sumItop - sumIbot));

        end
    end
    rough_scan_table = imgaussfilt(rough_scan_table,1);
    RST_weights = 1 - mat2gray(rough_scan_table_tbd);
    rough_scan_table_weighted = imgaussfilt((rough_scan_table .* RST_weights),1);
    [K_scan_center_i,E_scan_center_i] = find(rough_scan_table_weighted==max(rough_scan_table_weighted(:)));
    
    %figure, 
    %subplot(1,3,1), imagesc(rough_scan_table), axis xy, 
    %subplot(1,3,2), imagesc(rough_scan_table_tbd), axis xy
    %subplot(1,3,3), imagesc(rough_scan_table_weighted), axis xy, hold on;
    %plot(E_scan_center_i,K_scan_center_i,'r+'), hold off;


    %%%%%%%%%%Now the actual ABEK scanning%%%%%%%%%
    
    %Setting the refined-E,K scan range based off rough scan results
    K_off_range = K_rough_range(K_scan_center_i) + K_ref_range;
    E_off_range = E_rough_range(E_scan_center_i) + E_ref_range;   
    
    
    EK_AB_MC_table = zeros(size(length(K_off_range),length(E_off_range)));
    EK_AB_A_table = zeros(size(length(K_off_range),length(E_off_range)));
    EK_AB_B_table = zeros(size(length(K_off_range),length(E_off_range)));
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
            end
                    
            %Create corr table for AvsB at each E,K point       
            AB_MC_table = zeros(size(length(A_range),length(B_range)));
            A_IT_i = 1;
            for A_i = 1:length(A_range)
                A = A_range(A_i);
                for B_i = 1:length(B_range)
                    B = B_range(B_i);          
                        
                    fass_window = fass([1:size(draw_box,1)]+K_off, [1:size(draw_box,2)]+E_off);
                    fass_window_norm = window_processor(fass_window);
                    
                    AB_MC_table(A_i,B_i) = sum(dot(fass_window_norm, ABBA_ITs(:,:,A_IT_i)));
                    A_IT_i = A_IT_i+1;
                end
            end 
            %Choose A,B that gave max corr at the current E,K point, store in table
            [A_it_i,B_it_i] = find(AB_MC_table==max(AB_MC_table(:)),1);
            A_it = A_range(A_it_i);
            B_it = B_range(B_it_i);

            EK_AB_MC_table(K_off_i,E_off_i) = max(AB_MC_table(:));
            EK_AB_A_table(K_off_i,E_off_i) = A_it;
            EK_AB_B_table(K_off_i,E_off_i) = B_it;      
        end    
    end
    
    %Find E,k with max corr value, then ID which A,B values were used there
    [K_off_it_i, E_off_it_i] = find(EK_AB_MC_table==max(EK_AB_MC_table(:)));
    K_off_it = K_off_range(K_off_it_i);
    E_off_it = E_off_range(E_off_it_i);       
   
    A_it_it = EK_AB_A_table(K_off_it_i,E_off_it_i);
    B_it_it = EK_AB_B_table(K_off_it_i,E_off_it_i);
    MC_it = EK_AB_MC_table(K_off_it_i,E_off_it_i);
    
    EK_AB_MCn_table = (EK_AB_MC_table-min(EK_AB_MC_table(:)))./(max(EK_AB_MC_table(:))-min(EK_AB_MC_table(:)));
    [max_a,max_b] = find(EK_AB_MCn_table>=MC_TH);
    corr_spread = length(max_a);
   
    %Draw the chosen A,B curves and the chosen fass window
    A_IT_i = find(A_range == A_it_it);
    B_IT_i = find(B_range == B_it_it);
    IT_IT = ABBA_ITs(:,:,(A_IT_i-1)*length(B_range)+B_IT_i);
    
    fass_window_it = fass([1:size(draw_box,1)]+K_off_it, [1:size(draw_box,2)]+E_off_it);
    fass_window_norm_it = window_processor(fass_window_it); 
    
    %%%%%%%%%%Few things to judge quality%%%%%%%%%%%%%%%
    
    %%%Combadge is the slice above dirac point where it should be
    %%%low-intensity (btw dirac point and bottom of conducting band)
    combadge_yp = A_it_it*(draw_x0) + B_it_it*((draw_x0).^2) + E_0;
    combadge_yn = -A_it_it*(draw_x0) + B_it_it*((draw_x0).^2) + E_0;

    combadge_itp = horzcat(combadge_yp,draw_x0+K_0);
    combadge_itn = horzcat(combadge_yn,draw_x0+K_0);

    combadge_itp_curt = combadge_itp(max(1,round(-A_it_it/(2*B_it_it))+K_0):end,:);
    combadge_itn_curt = combadge_itn(1:min(length(draw_x0),round(K_0+A_it_it/(2*B_it_it))),:);

    combadge_itp_it = reshape(combadge_itp_curt',1,[]);
    combadge_itn_it = reshape(combadge_itn_curt',1,[]);

    combadge_ITP = insertShape(draw_box, 'Line', combadge_itp_it);
    combadge_ITN = insertShape(draw_box, 'Line', combadge_itn_it);
    combadge_ITP = combadge_ITP(:,:,1);
    combadge_ITN = combadge_ITN(:,:,1);
    
    combadge_IT = combadge_ITP + combadge_ITN; 
    combadge_IT = imgaussfilt(combadge_IT, krillin_sigma); %to not penalize intensity from a good surface band that's slightly broad
    combadge = combadge_IT;
    combadge(combadge~=0)=1;
    combadge = abs(1-combadge);

    starfleet = bwconncomp(combadge);
    combadge(starfleet.PixelIdxList{1})=0;
    combadge(starfleet.PixelIdxList{2})=0;
    combadge(starfleet.PixelIdxList{3})=0;
    
    combadge(:,E_0+25:end) = zeros(size(combadge,1), size(combadge,2) - E_0 - 24);
    
    incombadge = sum(dot(combadge,fass_window_norm_it))/nnz(combadge) * 1000; 
    
    %Difference in sum of intensities between k>0 half and k<0 half for
    %each E column of chosen fass window, averaged
    sumIbot = sum(fass_window_norm_it(1:K_0,1:E_0+10));
    sumItop = sum(fass_window_norm_it(K_0+1:end,1:E_0+10));
    sumIdiff = mean(abs(sumItop-sumIbot));

    
    %Create the AvsB scan table at the chosen E,k point for analysis
    ABatDP = zeros(length(A_range), length(B_range));
    ABatDP_i = 1;
    for aa_i = 1:length(A_range)
        for bb_i = 1:length(B_range)
            ABatDP(aa_i,bb_i) = sum(dot(fass_window_it,ABBA_ITs(:,:,ABatDP_i)));
            ABatDP_i = ABatDP_i+1;
        end
    end
    ABatDP_gray = (ABatDP - min(ABatDP(:))) / (max(ABatDP(:)) - min(ABatDP(:)));
    
            
    
    ABEK_Es(i) = (E_off_it + E_0) * bin_E + round(bin_E/2);
    %ABEK_ks(i) = (K_off_it + K_0 + fass_k_off) * bin_k + round(bin_k/2);
    ABEK_ks(i) = (K_off_it + K_0)*bin_k + round(bin_k/2);
    ABEK_As(i) = A_it_it;
    ABEK_Bs(i) = B_it_it;
    ABEK_MCs(i) = MC_it;
    ABEK_MCSs(i) = corr_spread;
    TBDs(i) = sumIdiff;
    combadges(i) = incombadge;
    
    %ABjudge(i) = sum(dot(ABatDP_gray, goodAB)) / (size(goodAB,1)*size(goodAB,2));  
    
    draw_pic_it = zeros(size(fass));
    draw_pic_it([1:size(draw_box,1)] + K_off_it, [1:size(draw_box,2)] + E_off_it ) = IT_IT;

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
            drawww_scan_window = fass([1:size(draw_box,1)]+K_conf, [1:size(draw_box,2)]+E_conf);
            drawww_scan_window_norm = window_processor(drawww_scan_window);

            smear_table(K_conf_i,E_conf_i) = sum(dot(drawww_scan_window_norm,ABBA_ITs(:,:,(A_IT_i-1)*length(B_range)+B_IT_i)));
        end
    end

    gray_smear_table = mat2gray(smear_table);
    the_smear = length(find(gray_smear_table >= multi_TH))/nnz(gray_smear_table);

    if wannasee == 1    
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
        %rough_scan_tablea = rot90(rough_scan_table,-1);
        %imagesc(rough_scan_tablea), axis xy
        imagesc(B_range,A_range,ABatDP), axis xy
        xlabel('square (B) term')
        ylabel('linear (A) term')
        colormap(ax2, jet)
        title('A vs B at DP','FontSize',8)
        
        ax1 = subplot(2,4,1);
        FixedWidth = get(0,'FixedWidthFontName');
        %text(0.2,0.5,{['A (eVA) :  ',num2str(A_it_it*bin_E/bin_k*.8107)];['B (eVA) :  ',num2str(B_it_it*bin_E/bin_k*.8107)];['MC:  ',num2str(MC_it)];...
        text(0,0.5,{...
            ['A(bpix):  ',num2str(A_it_it)];...
            ['B(bpix):  ',num2str(B_it_it)];...
            ['E(bpix):  ',num2str(E_off_it+E_0)];...
            ['K(bpix):  ',num2str(K_off_it+K_0)];...
            ['  '];...
            ['TBD: ',num2str(TBDs(i))];...
            ['combadge: ',num2str(combadges(i))]}, 'FontName',FixedWidth,'FontSize',8);
        axis off

        ax3 = subplot(2,4,6); 
        gray_smear_tablea = rot90(gray_smear_table,-1);
        imagesc(K_conff_range,E_conff_range,gray_smear_tablea), axis xy
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
        title(num2str(i));
        %colormap(ax6, winter)
        
        suptitle({['Scan i=',num2str(i)]})
        pause(.01)
    end
end

toc;



