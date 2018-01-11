box_eV_range = [+.45,+.40];%
box_invA_range = [-.30, -.25];%[.07,.12];%
be2 = round(box_eV_range(2)/pix2eV);
be1 = round(box_eV_range(1)/pix2eV);
bk1 = box_invA_range(1)/pix2invA;
bk2 = box_invA_range(2)/pix2invA;

bgnd_box_eV_range = [+0.6, +0.25];%
bgnd_box_invA_range = [-0.35, -0.20]; 
bgnd_be2 = bgnd_box_eV_range(2)/pix2eV;
bgnd_be1 = bgnd_box_eV_range(1)/pix2eV;
bgnd_bk1 = bgnd_box_invA_range(1)/pix2invA;
bgnd_bk2 = bgnd_box_invA_range(2)/pix2invA;

bgnd_Iperpix = zeros(1,961);
 for i = find(kLOS>0)
        bgnd_res1i = result1i(:,:,i);
        if kLOS(i) + bgnd_bk1 < 1
            bgnd_res1i_box = bgnd_res1i(1:round(kLOS(i)+bgnd_bk2),round(AVE_FL-bgnd_be1):round(AVE_FL-bgnd_be2));
        elseif kLOS(i) + bgnd_bk2 > size(bgnd_res1i,1)
            bgnd_res1i_box = bgnd_res1i(round(kLOS(i)+bgnd_bk1):end,round(AVE_FL-bgnd_be1):round(AVE_FL-bgnd_be2));
        else
            bgnd_res1i_box = bgnd_res1i(round(kLOS(i)+bgnd_bk1):round(kLOS(i)+bgnd_bk2), round(AVE_FL-bgnd_be1):round(AVE_FL-bgnd_be2));%res1i(round(kLOS(i)+bk1):round(kLOS(i)+bk2), round(AVE_FL-be1):round(AVE_FL-be2));
        end
        bgnd_number_of_pixels = size(bgnd_res1i_box,1)*size(bgnd_res1i_box,2);
        bgnd_Iperpix(i) =  sum(sum(bgnd_res1i_box))/bgnd_number_of_pixels;
end
ave_bgnd_Iperpix = mean(bgnd_Iperpix)


figure,
imagesc(imgaussfilt(result1i(:,:,random_scan),5)), axis xy, hold on;
if size(bgnd_res1i_box,1) == size(result1i(:,:,1),1)
    plot([round(AVE_FL-bgnd_be1),round(AVE_FL-bgnd_be1),round(AVE_FL-bgnd_be2),round(AVE_FL-bgnd_be2),round(AVE_FL-bgnd_be1)],...
    [1,size(res1i,1),size(res1i,1),1,1],'r', 'DisplayName','ave bgnd'), hold on;
else
    plot([round(AVE_FL-bgnd_be1),round(AVE_FL-bgnd_be1),round(AVE_FL-bgnd_be2),round(AVE_FL-bgnd_be2),round(AVE_FL-bgnd_be1)],...
    [round(kLOS(i)+bgnd_bk1),round(kLOS(i)+bgnd_bk2),round(kLOS(i)+bgnd_bk2),round(kLOS(i)+bgnd_bk1),round(kLOS(i)+bgnd_bk1)],'r','DisplayName','ave bgnd'), hold on;
end
title({['average Iperpix all scans = ',num2str(ave_bgnd_Iperpix)]; ['E range [',num2str(bgnd_box_eV_range(1)),', ',num2str(bgnd_box_eV_range(2)),']eV']});

ind_panel_MBI = zeros(1,size(region_list,1));
all_in_panels = zeros(1,961) - 100;
all_in_panels_no_subtract = zeros(1,961) - 100;
panel_nnn = zeros(1,size(region_list,1));
NPI = 0;

scan_box_events = zeros(1,961,size(region_list,1)) - 100;
for panel = 1:size(region_list,1)
    npi = 0;
    for p_i = find(involved_scans_panel(:,:,panel)>0)   
        npi = npi+1;
        res1i = (result1i(:,:,p_i));
        res1i_box = res1i(round(kLOS(p_i)+bk1):round(kLOS(p_i)+bk2),round(AVE_FL-be1):round(AVE_FL-be2));
        number_of_pixels = size(res1i_box,1)*size(res1i_box,2);
        events_in_box_total = sum(sum(res1i_box));
        events_in_box_minus_bgnd = events_in_box_total - number_of_pixels*bgnd_Iperpix(p_i);       
        scan_box_events(1,p_i,panel) = events_in_box_minus_bgnd;
        all_in_panels_no_subtract(p_i) = events_in_box_total;
        all_in_panels(p_i) = events_in_box_minus_bgnd;
    end
    current_panel = scan_box_events(:,:,panel);
    ind_panel_MBI(panel) = mean(current_panel(current_panel>-100));
    panel_nnn(panel) = npi;
    NPI = NPI + npi;
end
all_panels_MBI = mean(all_in_panels(all_in_panels>-100))

figure, 
ax = subplot(3,1,1);
random_scan = round(961*rand);
imagesc(imgaussfilt(result1i(:,:,random_scan),5)), axis xy, hold on;
if size(bgnd_res1i_box,1) == size(result1i(:,:,1),1)
    plot([round(AVE_FL-bgnd_be1),round(AVE_FL-bgnd_be1),round(AVE_FL-bgnd_be2),round(AVE_FL-bgnd_be2),round(AVE_FL-bgnd_be1)],...
    [1,size(res1i,1),size(res1i,1),1,1],'r', 'DisplayName','ave bgnd'), hold on;
else
    plot([round(AVE_FL-bgnd_be1),round(AVE_FL-bgnd_be1),round(AVE_FL-bgnd_be2),round(AVE_FL-bgnd_be2),round(AVE_FL-bgnd_be1)],...
    [round(kLOS(i)+bgnd_bk1),round(kLOS(i)+bgnd_bk2),round(kLOS(i)+bgnd_bk2),round(kLOS(i)+bgnd_bk1),round(kLOS(i)+bgnd_bk1)],'r','DisplayName','ave bgnd'), hold on;
end
plot([round(AVE_FL-be1),round(AVE_FL-be1),round(AVE_FL-be2),round(AVE_FL-be2),round(AVE_FL-be1)],...
    [round(kLOS(i)+bk1),round(kLOS(i)+bk2),round(kLOS(i)+bk2),round(kLOS(i)+bk1),round(kLOS(i)+bk1)],'w','DisplayName','set box'), hold off;
title({['Red for bgnd: [',num2str(bgnd_box_eV_range(1)),', ',num2str(bgnd_box_eV_range(2)),'] eV'];['White for box: [',num2str(box_eV_range(1)),', ',num2str(box_eV_range(2)),'] eV,   [',num2str(box_invA_range(1)),', ',num2str(box_invA_range(2)),'] invA']});
%title(['above-FL I per pix = ',num2str(ave_above_Iperpix),'  scan ',num2str(random_scan)])

ax = subplot(3,1,2);
h = histogram(all_in_panels_no_subtract(all_in_panels_no_subtract>-100),'BinWidth',1); 
xlabel('Events in white box')
ylabel('Scans')
title('Without subtracting background')

ax = subplot(3,1,3);
hh = histogram(all_in_panels(all_in_panels>-100),'BinWidth',1);
xlabel('Events in white box')
ylabel('Scans')
title('With subtracting background')
% for panel = 1:size(region_list,1)%:size(region_list,1)
%     panel_scans = scan_box_events(:,:,panel);
%     panel_scan_Is = panel_scans(panel_scans>-1);
% 
%     panel_scan_Ps = zeros(size(panel_scan_Is));
%     panel_scan_Ptots = zeros(size(panel_scan_Is));
%     if panel == 3
% %        disp('panel_scans='),disp(num2str(panel_scans))
% %         disp('panel_scan_Is='),disp(num2str(panel_scan_Is))
% %         disp('panel_scan_Ps='),disp(num2str(panel_scan_Ps))
% %         disp('panel_scan_Ptots='),disp(num2str(panel_scan_Ptots))
%     end
%     for pani = 1:length(panel_scan_Is)
%         panel_scan_Ps(pani) = (MBI ^ panel_scan_Is(pani)) * exp(-MBI) / factorial(panel_scan_Is(pani));
%         panel_scan_Ptots(pani) = (all_MBI ^ panel_scan_Is(pani)) * exp(-all_MBI) / factorial(panel_scan_Is(pani));
%     end
%     PBI_ind_ave(panel) = mean(panel_scan_Ps);
%     PBI_ind_ave_all(panel) = mean(panel_scan_Ptots);
% end
% 
% figure,
% II_n = 1;
% for II = 1:size(region_list,1)
%     spec = out_spec{1,II};
%     eax = out_eax{1,II};
%     kax = out_kax{1,II};
% 
%     PBI = round(PBIs(II));
%     nnnn = panel_nnn(II);
% 
%     poi = ((MBI^PBI)*(exp(-MBI))) / factorial(PBI);
%     poii = PBI_ind_ave(II);
%     poiii = PBI_ind_ave_all(II);
% 
%     ax = subplot(2,size(region_list,1),II_n);
%     imagesc(kax, flip(eax), rot90(norman(spec,0,5))), axis xy, hold on;
%     plot([box_invA_range(1),box_invA_range(2),box_invA_range(2),box_invA_range(1),box_invA_range(1)],-[boxtop_eV,boxtop_eV,boxbot_eV,boxbot_eV,boxtop_eV],'r'), hold off;
%     descr = {['BCB:',num2str(B_interval_list_abs(II,1)),'-',num2str(B_interval_list_abs(II,2))]};
%     axes(ax)
%     text(-.15,-.85,descr,'FontSize',8)    
%     title({['PBI=',num2str(PBI)];['P\_ave=',num2str(round(poiii,4))];['nnn=',num2str(nnnn)]},'FontSize',8)
%     ylim([-.7,0.25])
%     II_n = II_n+1;
%     colormap jet
% end
% 
% ax = subplot(2,size(region_list,1),[size(region_list,1)+1 2*size(region_list,1)]);
% yyaxis right 
% 
% h = histogram(all_panels(all_panels>-1),BI_edges); 
% hold on;
% BI_binmids = h.BinEdges(1:end-1)+h.BinWidth/2;%h.BinEdges(1:BI_nbins)+h.BinWidth/2;
% BI_binvals = smooth(h.Values); %./ max(smooth(h.Values));
% BI_poissrat = zeros(size(BI_binvals));
% poi_x = [1:1:round(BI_binmids(end))];
% poi_yfactorials = zeros(size(poi_x));
% gaulim = zeros(size(poi_x));
% for ii = poi_x
%     poi_yfactorials(ii) = factorial(ii);
%     gaulim(ii) = exp((-(ii - MBI)^2)/(2*MBI)) / (sqrt(2*pi*MBI));
% end
% poi_y = ((MBI.^poi_x)*exp(-MBI)) ./ poi_yfactorials;
% %poi_normy = poi_y / max(poi_y);
% if max(poi_y) == 0
%     disp('Max poisson value is 0')
% elseif max(poi_y) == inf
%     disp('Max poisson value is INF')
% end
% BI_histarea = (h.BinWidth * sum(h.Values));
% area_norm = BI_histarea / trapz(poi_y);
% for BIi = 1:(length(BI_edges)-1)
%     if ceil(BI_binmids(BIi))> max(poi_x)
%         continue
%     end
%     BI_val = BI_binvals(BIi);
%     if BI_val == 0
%         BI_poissrat(BIi) = 0;
%         continue
%     end
%     poissval_compare = poi_y(find(poi_x==ceil(BI_binmids(BIi))));
%     %BI_poissrat(BIi) = BI_val / poissval_compare;
%     BI_poissrat(BIi) = area_norm*poissval_compare / BI_val;
% end  
% 
% plot(poi_x, (BI_histarea / trapz(gaulim))*gaulim,'b', 'DisplayName','Gau lim'), hold on;
% plot(poi_x,poi_y,'DisplayName','Pois abt (inv)MBI'), hold off;
% yyaxis left 
% semilogy(BI_binmids,BI_poissrat,'-m','LineWidth',2, 'DisplayName', 'Ratio poi:hist');
% grid on
% legend('show')
% xlabel('box events "BI"','FontSize',9)
% title(['MBI (inv pans only) = ',num2str(round(MBI,1)),'    Box E: ',num2str(boxtop_eV),' - ',num2str(boxbot_eV),' eV     Box K: ',num2str(box_invA_range(1)),'-',num2str(box_invA_range(2)),' invA'],'FontSize',9)
% 
% conv_bgnd_poi = conv(ave_bgnd_dist_interped,poi_y);
% figure, 
% subplot(1,2,1)
% yyaxis left
% plot(poi_x,ave_bgnd_dist_interped,'r','DisplayName','ave bgnd above FL')
% yyaxis right
% plot(poi_x,poi_y,'b', 'DisplayName', 'Gauss abt MBI of current box')
% legend('show')
% 
% subplot(1,2,2),
% plot(conv_bgnd_poi,'m', 'DisplayName', 'Conv of avebgnd and Gauss abt MBI'), hold on,
% plot(BI_binmids,BI_binvals,'c'), hold off

