% 
% box_eV_ranges =[-.2,-.15;   -.175,-.125;    -.15,-.1; -.125, -.075];% [.45,.5];%
% box_invA_ranges = [-.12,-.07;   -.025,.025;   .07,.12];% [.08,.13];%
% 
% box_N = 0;
% BI_edges = 0:2:40;%
% ave_bgnd_dist = zeros(size(1,length(BI_edges)-1));
% for box_eV_range_i = 1:size(box_eV_ranges,1)
%     box_eV_range = box_eV_ranges(box_eV_range_i,:);
%     for box_invA_range_i = 1:size(box_invA_ranges,1)
%         box_invA_range = box_invA_ranges(box_invA_range_i,:);
%         box_N = box_N + 1;
% 
%         boxtop_eV = box_eV_range(1);%.24;% 0.2;
%         boxbot_eV = box_eV_range(2);%.30;%0.30;
% 
%         be2 = round(boxtop_eV/pix2eV);
%         be1 = round(boxbot_eV/pix2eV);
% 
%         bk1 = box_invA_range(1)/pix2invA;
%         bk2 = box_invA_range(2)/pix2invA;
% 
% 
%         PBIs = zeros(1,size(region_list,1));
%         PBI_ind_ave = zeros(size(PBIs));
%         panel_nnn = zeros(1,size(region_list,1));
% 
%         all_panels = zeros(1,961) - 100;
%         all_hi_dpi = zeros(1,961);
% 
% 
%          for i = find(kLOS>0)
%                 res1i = result1i(:,:,i);
%                 all_hi_dpi(i) = round(sum(sum(res1i(round(kLOS(i)+bk1):round(kLOS(i)+bk2), round(AVE_FL-be1):round(AVE_FL-be2)))));%round(rfc_FL_Es(i)-be1):round(rfc_FL_Es(i)-be2)))));
%         end
%         all_MBI = mean(all_hi_dpi(all_hi_dpi>0));
%         panel_mean = zeros(1,size(region_list,1));
%         all_panels_mean = 0;
%         NPI = 0;
%         scan_box_events = zeros(1,961,size(region_list,1)) - 100;
%         for panel = 1:size(region_list,1)
%             %scan_box_events = zeros(size(find(involved_scans_panel(:,:,panel)>0)));
%             %scan_box_events = zeros(1,961);
%             %disp('scan_box_events size'),disp(num2str(size(scan_box_events)))
%             %scan_box_probs = zeros(size(find(involved_scans_panel(:,:,panel)>0)));
%             npi = 0;
%             for p_i = find(involved_scans_panel(:,:,panel)>0)   
%                 npi = npi+1;
%                 res1i = (result1i(:,:,p_i));
%                 current_events = round(sum(sum(res1i(round(kLOS(p_i)+bk1):round(kLOS(p_i)+bk2), round(AVE_FL-be1):round(AVE_FL-be2)))));%round(rfc_FL_Es(i)-be1):round(rfc_FL_Es(i)-be2)))));
%                % scan_box_events(p_i,panel) = current_events;
%                 scan_box_events(1,p_i,panel) = current_events;
% 
%                 panel_mean(panel) = panel_mean(panel) + current_events;
%                 all_panels_mean = all_panels_mean + current_events;
%                 all_panels(p_i) = current_events;
%             end
%             %disp(num2str(npi))
%             %panel_scan_Is = scan_box_events(scan_box_events>0);
%             %panel_scan_Is = panel_scans(panel_scans>0);
%             PBIs(panel) = (panel_mean(panel)/npi);
%             %PBIs(panel) = mean(scan_box_events(scan_box_events>0))
% 
% 
%             panel_nnn(panel) = npi;
%             NPI = NPI + npi;%length(panel_scan_Is)
%             %panel_scan_Ps = zeros(size(panel_scan_Is));
%             %panel_scan_Ptots = zeros(size(panel_scan_Is));
% 
%             %figure, histogram(panel_scan_Is,20), title(['Panel#',num2str(panel)])
%         %     for pani = 1:length(panel_scan_Is)
%         %         panel_scan_Ps(pani) = (PBIs(panel) ^ panel_scan_Is(pani)) * exp(-PBIs(panel)) / factorial(panel_scan_Is(pani));
%         %     end
%         %     PBI_ind_ave(panel) = mean(panel_scan_Ps);
%         %     figure, histogram(panel_scan_Ps)
%         end
%         panels_mean = mean(PBIs);
%         MBI = all_panels_mean / NPI
%         %MBI = mean(all_panels(all_panels>0))%panels_mean;
% 
% 
%         for panel = 1:size(region_list,1)%:size(region_list,1)
%             panel_scans = scan_box_events(:,:,panel);
%             panel_scan_Is = panel_scans(panel_scans>-1);
% 
%             panel_scan_Ps = zeros(size(panel_scan_Is));
%             panel_scan_Ptots = zeros(size(panel_scan_Is));
%             if panel == 3
%         %        disp('panel_scans='),disp(num2str(panel_scans))
%         %         disp('panel_scan_Is='),disp(num2str(panel_scan_Is))
%         %         disp('panel_scan_Ps='),disp(num2str(panel_scan_Ps))
%         %         disp('panel_scan_Ptots='),disp(num2str(panel_scan_Ptots))
% 
%             end
%             for pani = 1:length(panel_scan_Is)
%                 panel_scan_Ps(pani) = (MBI ^ panel_scan_Is(pani)) * exp(-MBI) / factorial(panel_scan_Is(pani));
%                 panel_scan_Ptots(pani) = (all_MBI ^ panel_scan_Is(pani)) * exp(-all_MBI) / factorial(panel_scan_Is(pani));
%             end
%             PBI_ind_ave(panel) = mean(panel_scan_Ps);
%             PBI_ind_ave_all(panel) = mean(panel_scan_Ptots);
%         end
% 
%         figure,
% 
%         II_n = 1;
%         for II = 1:size(region_list,1)
%             spec = out_spec{1,II};
%             eax = out_eax{1,II};
%             kax = out_kax{1,II};
% 
%             PBI = round(PBIs(II));
%             nnnn = panel_nnn(II);
% 
%             poi = ((MBI^PBI)*(exp(-MBI))) / factorial(PBI);
%             poii = PBI_ind_ave(II);
%             poiii = PBI_ind_ave_all(II);
% 
%             ax = subplot(2,size(region_list,1),II_n);
%             imagesc(kax, flip(eax), rot90(norman(spec,0,5))), axis xy, hold on;
%             plot([box_invA_range(1),box_invA_range(2),box_invA_range(2),box_invA_range(1),box_invA_range(1)],-[boxtop_eV,boxtop_eV,boxbot_eV,boxbot_eV,boxtop_eV],'r'), hold off;
%             descr = {['BCB:',num2str(B_interval_list_abs(II,1)),'-',num2str(B_interval_list_abs(II,2))]};
%             axes(ax)
%             text(-.15,-.85,descr,'FontSize',8)    
%             title({['PBI=',num2str(PBI)];['P\_ave=',num2str(round(poiii,4))];['nnn=',num2str(nnnn)]},'FontSize',8)
%             ylim([-.7,0.25])
%             II_n = II_n+1;
%             colormap jet
%         end
%         %suptitle({['MBI=',num2str(MBI),' boxE: ',num2str(boxtop_eV),'-',num2str(boxbot_eV),'  boxK: +/-',num2str(box_khalfwidth_invA)]})
% 
%         %all_MBI = 30;
% 
%         ax = subplot(2,size(region_list,1),[size(region_list,1)+1 2*size(region_list,1)]);
%         yyaxis right 
%         %BI_nbins = 20;
%         
%         h = histogram(all_panels(all_panels>-1),BI_edges); 
%         hold on;
%         BI_binmids = h.BinEdges(1:end-1)+h.BinWidth/2;%h.BinEdges(1:BI_nbins)+h.BinWidth/2;
%         BI_binvals = smooth(h.Values); %./ max(smooth(h.Values));
%         BI_poissrat = zeros(size(BI_binvals));
%         poi_x = [1:1:round(BI_binmids(end))];
%         poi_yfactorials = zeros(size(poi_x));
%         gaulim = zeros(size(poi_x));
%         for ii = poi_x
%             poi_yfactorials(ii) = factorial(ii);
%             gaulim(ii) = exp((-(ii - all_MBI)^2)/(2*all_MBI)) / (sqrt(2*pi*all_MBI));
%         end
%         poi_y = ((all_MBI.^poi_x)*exp(-all_MBI)) ./ poi_yfactorials;
%         %poi_normy = poi_y / max(poi_y);
%         if max(poi_y) == 0
%             disp('Max poisson value is 0')
%         elseif max(poi_y) == inf
%             disp('Max poisson value is INF')
%         end
%         BI_histarea = (h.BinWidth * sum(h.Values));
%         area_norm = BI_histarea / trapz(poi_y);
%         for BIi = 1:(length(BI_edges)-1)
%             if ceil(BI_binmids(BIi))> max(poi_x)
%                 continue
%             end
%             BI_val = BI_binvals(BIi);
%             if BI_val == 0
%                 BI_poissrat(BIi) = 0;
%                 continue
%             end
%             poissval_compare = poi_y(find(poi_x==ceil(BI_binmids(BIi))));
%             %BI_poissrat(BIi) = BI_val / poissval_compare;
%             BI_poissrat(BIi) = area_norm*poissval_compare / BI_val;
%         end  
% 
% 
%         %plot(poi_x, area_norm*poi_y, 'r'), legend('Poisson'), 
%         %hold on;
%         plot(poi_x, (BI_histarea / trapz(gaulim))*gaulim,'b'), legend('Gaussian'), hold off;
%         %plot([0,poi_x(end)],[1,1],'--m'), hold off;
% 
%         yyaxis left 
%         semilogy(BI_binmids,BI_poissrat,'-m','LineWidth',2);
%         grid on
%         legend('Ratio poi:bin')
%         legend('show')
%         %ylim([0 BI_poissrat(
%         xlabel('box events "BI"','FontSize',9)
%         title(['allMBI = ',num2str(round(all_MBI,1)),'(allinvMBI=',num2str(round(MBI,1)),')     Box E: ',num2str(boxtop_eV),' - ',num2str(boxbot_eV),' eV     Box K: ',num2str(box_invA_range(1)),'-',num2str(box_invA_range(2)),' invA'],'FontSize',9)
%         %title({['MBI = ',num2str(MBI),' events'];['Box E range: ',num2str(boxtop_eV),' - ',num2str(boxbot_eV),' eV'];['Box K range: +/- ',num2str(box_khalfwidth_invA),' invA']})
% 
%         %{
%         for panel = 1:size(region_list,1)%:size(region_list,1)
%             panel_scans = scan_box_events(:,:,panel);
%             panel_scan_Is = panel_scans(panel_scans>0);
%             figure, yyaxis left, histogram(panel_scan_Is,20), hold on;
%             yyaxis right, plot(poi_x,poi_normy), hold off;
%             title(['Panel#',num2str(panel),' for E range ',num2str(box_eV_range),' eV'])
%         end
%         %}   
%         ave_bgnd_dist = ave_bgnd_dist + BI_binvals;
%         
%     end
% end
% ave_bgnd_dist = ave_bgnd_dist ./ box_N;
% figure, plot(BI_binmids, ave_bgnd_dist,'r')
% 
% ave_bgnd_dist_interped = interp1(BI_binmids, ave_bgnd_dist, poi_x);



box_eV_range = [-.150,-.100];%
box_invA_range = [.07,.12];%

boxtop_eV = box_eV_range(1);%.24;% 0.2;
boxbot_eV = box_eV_range(2);%.30;%0.30;

be2 = round(boxtop_eV/pix2eV);
be1 = round(boxbot_eV/pix2eV);

bk1 = box_invA_range(1)/pix2invA;
bk2 = box_invA_range(2)/pix2invA;

PBIs = zeros(1,size(region_list,1));
PBI_ind_ave = zeros(size(PBIs));
panel_nnn = zeros(1,size(region_list,1));

all_panels = zeros(1,961) - 100;
all_panels_minus_bgnd = zeros(1,961) - 100;
all_hi_dpi = zeros(1,961);

 for i = find(kLOS>0)
        res1i = result1i(:,:,i);
        if ((kLOS(i) + bk1) < 0) || ((kLOS(i)+bk2) > 300)
            continue
        end
        all_hi_dpi(i) = round(sum(sum(res1i(round(kLOS(i)+bk1):round(kLOS(i)+bk2), ...
                                    round(AVE_FL-be1):round(AVE_FL-be2)))));
end
all_MBI = mean(all_hi_dpi(all_hi_dpi>0));
panel_mean = zeros(1,size(region_list,1));
all_panels_mean = 0;
NPI = 0;

scan_box_events = zeros(1,961,size(region_list,1)) - 100;
for panel = 1:size(region_list,1)
    npi = 0;
    for p_i = find(involved_scans_panel(:,:,panel)>0)   
        npi = npi+1;
        res1i = (result1i(:,:,p_i));
        box_pixnum = size(res1i,1)*size(res1i,2);
        current_events = round(sum(sum(res1i(round(kLOS(p_i)+bk1):round(kLOS(p_i)+bk2), round(AVE_FL-be1):round(AVE_FL-be2)))));%round(rfc_FL_Es(i)-be1):round(rfc_FL_Es(i)-be2)))));
       % scan_box_events(p_i,panel) = current_events;
        scan_box_events(1,p_i,panel) = current_events;

        panel_mean(panel) = panel_mean(panel) + current_events;
        all_panels_mean = all_panels_mean + current_events;
        all_panels(p_i) = current_events;
        all_panels_minus_bgnd(p_i) = current_events - box_pixnum*hypoFL_aveIperpix;
    end
    %disp(num2str(npi))
    %panel_scan_Is = scan_box_events(scan_box_events>0);
    %panel_scan_Is = panel_scans(panel_scans>0);
    PBIs(panel) = (panel_mean(panel)/npi);
    %PBIs(panel) = mean(scan_box_events(scan_box_events>0))


    panel_nnn(panel) = npi;
    NPI = NPI + npi;%length(panel_scan_Is)
end
panels_mean = mean(PBIs);
MBI = all_panels_mean / NPI
%MBI = mean(all_panels(all_panels>0))%panels_mean;

for panel = 1:size(region_list,1)%:size(region_list,1)
    panel_scans = scan_box_events(:,:,panel);
    panel_scan_Is = panel_scans(panel_scans>-1);

    panel_scan_Ps = zeros(size(panel_scan_Is));
    panel_scan_Ptots = zeros(size(panel_scan_Is));
    if panel == 3
%        disp('panel_scans='),disp(num2str(panel_scans))
%         disp('panel_scan_Is='),disp(num2str(panel_scan_Is))
%         disp('panel_scan_Ps='),disp(num2str(panel_scan_Ps))
%         disp('panel_scan_Ptots='),disp(num2str(panel_scan_Ptots))
    end
    for pani = 1:length(panel_scan_Is)
        panel_scan_Ps(pani) = (MBI ^ panel_scan_Is(pani)) * exp(-MBI) / factorial(panel_scan_Is(pani));
        panel_scan_Ptots(pani) = (all_MBI ^ panel_scan_Is(pani)) * exp(-all_MBI) / factorial(panel_scan_Is(pani));
    end
    PBI_ind_ave(panel) = mean(panel_scan_Ps);
    PBI_ind_ave_all(panel) = mean(panel_scan_Ptots);
end

figure,
II_n = 1;
for II = 1:size(region_list,1)
    spec = out_spec{1,II};
    eax = out_eax{1,II};
    kax = out_kax{1,II};

    PBI = round(PBIs(II));
    nnnn = panel_nnn(II);

    poi = ((MBI^PBI)*(exp(-MBI))) / factorial(PBI);
    poii = PBI_ind_ave(II);
    poiii = PBI_ind_ave_all(II);

    ax = subplot(2,size(region_list,1),II_n);
    imagesc(kax, flip(eax), rot90(norman(spec,0,5))), axis xy, hold on;
    plot([box_invA_range(1),box_invA_range(2),box_invA_range(2),box_invA_range(1),box_invA_range(1)],-[boxtop_eV,boxtop_eV,boxbot_eV,boxbot_eV,boxtop_eV],'r'), hold off;
    descr = {['BCB:',num2str(B_interval_list_abs(II,1)),'-',num2str(B_interval_list_abs(II,2))]};
    axes(ax)
    text(-.15,-.85,descr,'FontSize',8)    
    title({['PBI=',num2str(PBI)];['P\_ave=',num2str(round(poiii,4))];['nnn=',num2str(nnnn)]},'FontSize',8)
    ylim([-.7,0.25])
    II_n = II_n+1;
    colormap jet
end

ax = subplot(2,size(region_list,1),[size(region_list,1)+1 2*size(region_list,1)]);
yyaxis right 

h = histogram(all_panels(all_panels>-1),BI_edges); 
hold on;
BI_binmids = h.BinEdges(1:end-1)+h.BinWidth/2;%h.BinEdges(1:BI_nbins)+h.BinWidth/2;
BI_binvals = smooth(h.Values); %./ max(smooth(h.Values));
BI_poissrat = zeros(size(BI_binvals));
poi_x = [1:1:round(BI_binmids(end))];
poi_yfactorials = zeros(size(poi_x));
gaulim = zeros(size(poi_x));
for ii = poi_x
    poi_yfactorials(ii) = factorial(ii);
    gaulim(ii) = exp((-(ii - MBI)^2)/(2*MBI)) / (sqrt(2*pi*MBI));
end
poi_y = ((MBI.^poi_x)*exp(-MBI)) ./ poi_yfactorials;
%poi_normy = poi_y / max(poi_y);
if max(poi_y) == 0
    disp('Max poisson value is 0')
elseif max(poi_y) == inf
    disp('Max poisson value is INF')
end
BI_histarea = (h.BinWidth * sum(h.Values));
area_norm = BI_histarea / trapz(poi_y);
for BIi = 1:(length(BI_edges)-1)
    if ceil(BI_binmids(BIi))> max(poi_x)
        continue
    end
    BI_val = BI_binvals(BIi);
    if BI_val == 0
        BI_poissrat(BIi) = 0;
        continue
    end
    poissval_compare = poi_y(find(poi_x==ceil(BI_binmids(BIi))));
    %BI_poissrat(BIi) = BI_val / poissval_compare;
    BI_poissrat(BIi) = area_norm*poissval_compare / BI_val;
end  

plot(poi_x, (BI_histarea / trapz(gaulim))*gaulim,'b', 'DisplayName','Gau lim'), hold on;
plot(poi_x,poi_y,'DisplayName','Pois abt (inv)MBI'), hold off;
yyaxis left 
semilogy(BI_binmids,BI_poissrat,'-m','LineWidth',2, 'DisplayName', 'Ratio poi:hist');
grid on
legend('show')
xlabel('box events "BI"','FontSize',9)
title(['MBI (inv pans only) = ',num2str(round(MBI,1)),'    Box E: ',num2str(boxtop_eV),' - ',num2str(boxbot_eV),' eV     Box K: ',num2str(box_invA_range(1)),'-',num2str(box_invA_range(2)),' invA'],'FontSize',9)

conv_bgnd_poi = conv(ave_bgnd_dist_interped,poi_y);
figure, 
subplot(1,2,1)
yyaxis left
plot(poi_x,ave_bgnd_dist_interped,'r','DisplayName','ave bgnd above FL')
yyaxis right
plot(poi_x,poi_y,'b', 'DisplayName', 'Gauss abt MBI of current box')
legend('show')

subplot(1,2,2),
plot(conv_bgnd_poi,'m', 'DisplayName', 'Conv of avebgnd and Gauss abt MBI'), hold on,
plot(BI_binmids,BI_binvals,'c'), hold off

