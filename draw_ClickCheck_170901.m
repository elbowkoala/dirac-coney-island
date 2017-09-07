map_to_click = reshape(draw_MCSs_906,31,31); %%Input your map here%%%
dirac_Es = draw_Es_906;   %%%%Input your E vector (units of frame pixels)
dirac_ks = draw_ks_906;   %%%%Input your k vector

%%%%%%%%%%%%%
click_map = figure;

abcd = 1;

while abcd==1
    figure(click_map), imagesc(map_to_click), colormap jet, title('Click on the scan to check'), axis xy;
    [click_x, click_y] = ginput(1);
    click_xr = round(click_x);
    click_yr = round(click_y);
    clicked_i = (click_xr-1)*31 + click_yr;
    
 
    cone = cones(:,:,clicked_i);

    [fass, fass_k_off] = kLOSfinder5(cones(:,:,clicked_i),bin_E,bin_k);
    fass = imgaussfilt(fass,fass_sigma);
    %%%%%%%%%%
    E_off_ABBAn_it = ((draw_Es_906(clicked_i) - round(bin_E/2))/bin_E) - E_0;
    K_off_ABBAn_it = ((draw_ks_906(clicked_i) - round(bin_k/2))/bin_k) - K_0 - fass_k_off;
    ABBAn_it_it = draw_As_906(clicked_i);
    BAABn_it_it = draw_Bs_906(clicked_i);
    MCS_ABBAn_it = draw_MCSs_906(clicked_i);
    %%%%%%%%%%%%
    draw_yp_it_it = ABBAn_it_it*(draw_x-K_0) + BAABn_it_it*((draw_x-K_0).^2) + E_0;
    draw_yn_it_it = -ABBAn_it_it*(draw_x-K_0) + BAABn_it_it*((draw_x-K_0).^2) + E_0;
    
    draw_ittp = horzcat(draw_yp_it_it,draw_x);
    draw_ittn = horzcat(draw_yn_it_it,draw_x);
    
    draw_ittp_curt = draw_ittp(max(1,round(-ABBAn_it_it/(2*BAABn_it_it)+K_0)):end,:);
    draw_ittn_curt = draw_ittn(1:min(length(draw_x),round(K_0+ABBAn_it_it/(2*BAABn_it_it))),:);
    
    draw_ittp_it = reshape(draw_ittp_curt',1,[]);
    draw_ittn_it = reshape(draw_ittn_curt',1,[]);
    
    ITTP = insertShape(draw_box, 'Line', draw_ittp_it, 'LineWidth',LineWidth);
    ITTN = insertShape(draw_box, 'Line', draw_ittn_it, 'LineWidth',LineWidth);

    ITTP = mat2gray(ITTP(:,:,1));
    ITTN = mat2gray(ITTN(:,:,1));

    ITTP = imgaussfilt(ITTP,draw_sigma);
    ITTN = imgaussfilt(ITTN,draw_sigma);

    ITT = ITTP + ITTN;

    ITT = (ITT - min(ITT(:))) ./ (max(ITT(:))-min(ITT(:)));
    %ITT = ITT ./ sum(ITT(:));
    ITT_cropped = ITT(1+draw_margin:end-draw_margin,1+draw_margin:end-draw_margin);

    draw_scan_window_it = fass([1+draw_margin:size(draw_box,1)-draw_margin]+K_off_ABBAn_it,...
                            [1+draw_margin:size(draw_box,2)-draw_margin]+E_off_ABBAn_it);
    draw_scan_window_norm_it = (draw_scan_window_it - min(draw_scan_window_it(:)))./(max(draw_scan_window_it(:))-min(draw_scan_window_it(:))); %draw_scan_window_it ./ sum(draw_scan_window_it(:)); 
    
    draw_scan_window_norm_it = draw_scan_window_norm_it - mean(draw_scan_window_norm_it(:));
    
    draw_pic_it = zeros(size(fass));
    draw_pic_it([1+draw_margin:size(draw_box,1)-draw_margin] + K_off_ABBAn_it,...
             [1+draw_margin:size(draw_box,2)-draw_margin] + E_off_ABBAn_it ) = ITT_cropped;
    %{
    ABBA_norm_table_it = zeros(size(length(A_range),length(B_range)));
    for A_i = 1:length(A_range)
        A = A_range(A_i);

        for B_i = 1:length(B_range)
            B = B_range(B_i);

            draw_picc = zeros(size(fass));

            draw_ypp = A*(draw_x-K_0) + B*((draw_x-K_0).^2) + E_0;
            draw_ynn = -A*(draw_x-K_0) + B*((draw_x-K_0).^2) + E_0;

            draw_itpp = horzcat(draw_ypp,draw_x);
            draw_itnn = horzcat(draw_ynn,draw_x);

            draw_itpp_curt = draw_itpp(max(1,round(-A/(2*B)+K_0)):end,:);
            draw_itnn_curt = draw_itnn(1:min(length(draw_x),round(K_0+A/(2*B))),:);

            draw_itpp_it = reshape(draw_itpp_curt',1,[]);
            draw_itnn_it = reshape(draw_itnn_curt',1,[]);
            %{
            draw_itpp_it = zeros(1,2*length((max(1,round(-A/(2*B)+K_0)):length(draw_x))));
            itpp_i = 1;
            for ii = max(1,round(-A/(2*B)+K_0)):length(draw_x)
                draw_itpp_it(2*itpp_i-1:2*itpp_i) = draw_itpp(ii,:);
                itpp_i = itpp_i + 1;
            end

            draw_itnn_it = zeros(1,2*length(1:min(length(draw_x),round(K_0+A/(2*B)))));
            itnn_i = 1;
            for ii = 1:min(length(draw_x),round(K_0+A/(2*B)))
                draw_itnn_it(2*itnn_i-1:2*itnn_i) = draw_itnn(ii,:);
                itnn_i = itnn_i + 1;
            end
            %}
            ITPP = insertShape(draw_box, 'Line', draw_itpp_it, 'LineWidth',LineWidth);
            ITNN = insertShape(draw_box, 'Line', draw_itnn_it, 'LineWidth',LineWidth);

            ITPP = mat2gray(ITPP(:,:,1));
            ITNN = mat2gray(ITNN(:,:,1));

            ITPP = imgaussfilt(ITPP,draw_sigma);
            ITNN = imgaussfilt(ITNN,draw_sigma);

            ITTT = ITPP + ITNN;
            %ITTT_a = max(ITTT);
            %ITTT_b = repmat(ITTT_a,size(ITTT,1),1);
            %ITTT = ITTT ./ ITTT_b;

            ITTT = (ITTT - min(ITTT(:)))./(max(ITTT(:))-min(ITTT(:)));
            %ITTT = ITTT ./ sum(ITTT(:));

            ITTT_cropped = ITTT(1+draw_margin:end-draw_margin,1+draw_margin:end-draw_margin);

            draw_scan_window_itt = fass([1+draw_margin:size(draw_box,1)-draw_margin]+K_off_ABBAn_it,...
                            [1+draw_margin:size(draw_box,2)-draw_margin]+E_off_ABBAn_it);
            draw_scan_window_norm_itt = (draw_scan_window_itt-min(draw_scan_window_itt(:)))./(max(draw_scan_window_itt(:))-min(draw_scan_window_itt(:))); %draw_scan_window_itt ./ sum(draw_scan_window_itt(:));  
            %%%%%%
            draw_scan_window_norm_itt = draw_scan_window_norm_itt - mean(draw_scan_window_norm_itt(:));
            %%%%%
            draw_picc([1+draw_margin:size(draw_box,1)-draw_margin] + K_off_ABBAn_it,...
                     [1+draw_margin:size(draw_box,2)-draw_margin] + E_off_ABBAn_it ) = ITTT_cropped;

            ABBA_norm_table_it(A_i,B_i) = sum(sum(draw_scan_window_norm_itt.*ITTT_cropped));

        end
    end    
    %}
    figure,
    
    ax1 = subplot(2,1,1);
    E = mat2gray(fass);
    I = mat2gray(draw_pic_it);
    greenlol = cat(3,ones(size(E)),zeros(size(E)),zeros(size(E)));
    imshow(E,'InitialMag','fit','DisplayRange',[0 .5]), hold on;
    h = imshow(greenlol); 
    set(h,'AlphaData',I);
    hold off;
    colormap(ax1, winter)
    
    ax2 = subplot(2,1,2);
    FixedWidth = get(0,'FixedWidthFontName');
    %text(0.2,0.5,{['A (eVA) :  ',num2str(A_it_it*bin_E/bin_k*.8107)];['B (eVA) :  ',num2str(B_it_it*bin_E/bin_k*.8107)];['MC:  ',num2str(MC_it)];...
    text(0.2,0.5,{...
        ['          i = ',num2str(clicked_i)];...
        ['A (bpix):         ',num2str(ABBAn_it_it)];...
        ['B (bpix):         ',num2str(BAABn_it_it)];...
        ['MCS :             ',num2str(MCS_ABBAn_it)];...
        ['DP E coor (bpix): ',num2str(E_off_ABBAn_it+E_0)];...
        ['DP K coor (bpix): ',num2str(K_off_ABBAn_it+K_0)]} ,'FontName',FixedWidth);
    axis off
    
    %{
    ax5 = subplot(3,2,5);
    imagesc(draw_scan_window_it), axis xy
    colormap(ax5, pink)
    title('Best Window')

    ax6 = subplot(3,2,6);
    imagesc(ITT_cropped), axis xy
    colormap(ax6, bone)
    title('Best Draw Fit')

    ax3 = subplot(3,2,3); 
    imagesc(EK_ABBAn_table), axis xy
    colormap(ax3, jet)
    title('Max Corr (fass window normed)')
    %yticks([1,length(K_off_range)])
    %yticklabels({num2str(K_off_range(1)+K_0),num2str(K_off_range(end)+K_0)})
    %ylabel('DP K (binned pixels, fass frame)','FontSize',8)
    %xticks([1,length(E_off_range)])
    %xticklabels({num2str(E_off_range(1)+E_0),num2str(E_off_range(end)+E_0)});
    %xlabel('DP E (binned pixels, fass frame)','FontSize',8)

    ax4 = subplot(3,2,4);
    imagesc(ABBA_norm_table_it), axis xy, title('MaxCorr for A-vs-B at DP')
    colormap(ax4, jet)
    %{
    xticks([1:2:length(B_range)])
    xticklabels(num2str(B_range(1:2:length(B_range))))
    %}
    %set(gca, 'XTick',[1:2:length(B_range)], 'XTickLabel',{num2str(B_range_eVA(1,1:2:length(B_range))')}, 'FontSize',6)
    %set(gca, 'YTick',[1:2:length(A_range)], 'YTickLabel',{num2str(A_range_eVA(1,1:2:length(A_range))')}, 'FontSize',6)
    %set(gca, 'XTick',[1:2:length(B_range)], 'XTickLabel',{num2str(B_range(1,1:2:length(B_range))')}, 'FontSize',6)
    %set(gca, 'YTick',[1:2:length(A_range)], 'YTickLabel',{num2str(A_range(1,1:2:length(A_range))')}, 'FontSize',6)
    %yticks([1:2:length(A_range)])
    %yticklabels(num2str(A_range(1:2:length(A_range)')))
    %xlabel('B term (pixel units)','FontSize',8)
    %ylabel('A term (pixel units)','FontSize',8)
    %}
    
         
         
         
         
         
    
  
end