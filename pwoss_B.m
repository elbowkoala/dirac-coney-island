
box_eV_range = [.2,.33]-.4;
boxtop_eV = box_eV_range(1);%.24;% 0.2;
boxbot_eV = box_eV_range(2);%.30;%0.30;
%box_khalfwidth_invA = 0.1;
be2 = round(boxtop_eV/pix2eV);
be1 = round(boxbot_eV/pix2eV);
%bk = round(box_khalfwidth_invA/pix2invA);
box_invA_range = [-0.08, 0.02];
bk1 = box_invA_range(1)/pix2invA;
bk2 = box_invA_range(2)/pix2invA;


PBIs = zeros(1,size(region_list,1));
PBI_ind_ave = zeros(size(PBIs));
panel_nnn = zeros(1,size(region_list,1));

all_panels = zeros(1,961) - 100;
all_hi_dpi = zeros(1,961);


 for i = find(kLOS>0)
        res1i = result1i(:,:,i);
        all_hi_dpi(i) = round(sum(sum(res1i(round(kLOS(i)+bk1):round(kLOS(i)+bk2), round(AVE_FL-be1):round(AVE_FL-be2)))));%round(rfc_FL_Es(i)-be1):round(rfc_FL_Es(i)-be2)))));
end
all_MBI = mean(all_hi_dpi(all_hi_dpi>0));
panel_mean = zeros(1,size(region_list,1));
all_panels_mean = 0;
NPI = 0;
scan_box_events = zeros(1,961,size(region_list,1)) - 100;
for panel = 1:size(region_list,1)
    %scan_box_events = zeros(size(find(involved_scans_panel(:,:,panel)>0)));
    %scan_box_events = zeros(1,961);
    %disp('scan_box_events size'),disp(num2str(size(scan_box_events)))
    %scan_box_probs = zeros(size(find(involved_scans_panel(:,:,panel)>0)));
    npi = 0;
    for p_i = find(involved_scans_panel(:,:,panel)>0)   
        npi = npi+1;
        res1i = (result1i(:,:,p_i));
        current_events = round(sum(sum(res1i(round(kLOS(p_i)+bk1):round(kLOS(p_i)+bk2), round(AVE_FL-be1):round(AVE_FL-be2)))))%round(rfc_FL_Es(i)-be1):round(rfc_FL_Es(i)-be2)))));
       % scan_box_events(p_i,panel) = current_events;
        scan_box_events(1,p_i,panel) = current_events;
        
        panel_mean(panel) = panel_mean(panel) + current_events;
        all_panels_mean = all_panels_mean + current_events;
        all_panels(p_i) = current_events;
    end
    %disp(num2str(npi))
    %panel_scan_Is = scan_box_events(scan_box_events>0);
    %panel_scan_Is = panel_scans(panel_scans>0);
    PBIs(panel) = (panel_mean(panel)/npi);
    %PBIs(panel) = mean(scan_box_events(scan_box_events>0))
    

    panel_nnn(panel) = npi
    NPI = NPI + npi;%length(panel_scan_Is)
    %panel_scan_Ps = zeros(size(panel_scan_Is));
    %panel_scan_Ptots = zeros(size(panel_scan_Is));
    
    %figure, histogram(panel_scan_Is,20), title(['Panel#',num2str(panel)])
%     for pani = 1:length(panel_scan_Is)
%         panel_scan_Ps(pani) = (PBIs(panel) ^ panel_scan_Is(pani)) * exp(-PBIs(panel)) / factorial(panel_scan_Is(pani));
%     end
%     PBI_ind_ave(panel) = mean(panel_scan_Ps);
%     figure, histogram(panel_scan_Ps)
end
panels_mean = mean(PBIs)
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
%suptitle({['MBI=',num2str(MBI),' boxE: ',num2str(boxtop_eV),'-',num2str(boxbot_eV),'  boxK: +/-',num2str(box_khalfwidth_invA)]})


ax = subplot(2,size(region_list,1),[size(region_list,1)+1 2*size(region_list,1)]);
yyaxis right 
BI_nbins = 25;
h = histogram(all_panels(all_panels>-1),BI_nbins); 
hold on;
BI_binmids = h.BinEdges(1:BI_nbins)+h.BinWidth/2;
BI_binvals = smooth(h.Values); %./ max(smooth(h.Values));
BI_poissrat = zeros(size(BI_binvals));
poi_x = [1:1:round(BI_binmids(end))];
poi_yfactorials = zeros(size(poi_x));
gaulim = zeros(size(poi_x));
for ii = poi_x
    poi_yfactorials(ii) = factorial(ii);
    gaulim(ii) = exp((-(ii - all_MBI)^2)/(2*all_MBI)) / (sqrt(2*pi*all_MBI));
end
poi_y = ((all_MBI.^poi_x)*exp(-all_MBI)) ./ poi_yfactorials;
%poi_normy = poi_y / max(poi_y);
if max(poi_y) == inf || 0
    disp('Max poisson value is 0 OR INFINITY')
end

for BIi = 1:BI_nbins
    BI_val = BI_binvals(BIi)
    poissval_compare = poi_y(find(poi_x==round(BI_binmids(BIi))));
    %BI_poissrat(BIi) = BI_val / poissval_compare;
    BI_poissrat(BIi) = poissval_compare / BI_val;
end  
 
BI_histarea = (h.BinWidth * sum(h.Values));
plot(poi_x, (BI_histarea / trapz(poi_y))*poi_y, 'r'), legend('Poisson'), 
hold on;
plot(poi_x, (BI_histarea / trapz(gaulim))*gaulim,'b'), legend('Gaussian'), hold off;
%plot([0,poi_x(end)],[1,1],'--m'), hold off;

yyaxis left 
semilogy(BI_binmids,BI_poissrat,'-m','LineWidth',2);
legend('Ratio poi:bin')
legend('show')
%ylim([0 BI_poissrat(
xlabel('box events "BI"','FontSize',9)
title(['allMBI = ',num2str(round(all_MBI,1)),'(allinvMBI=',num2str(round(MBI,1)),')     Box E: ',num2str(boxtop_eV),' - ',num2str(boxbot_eV),' eV     Box K: ',num2str(box_invA_range(1)),'-',num2str(box_invA_range(2)),' invA'],'FontSize',9)
%title({['MBI = ',num2str(MBI),' events'];['Box E range: ',num2str(boxtop_eV),' - ',num2str(boxbot_eV),' eV'];['Box K range: +/- ',num2str(box_khalfwidth_invA),' invA']})

%{
for panel = 1:size(region_list,1)%:size(region_list,1)
    panel_scans = scan_box_events(:,:,panel);
    panel_scan_Is = panel_scans(panel_scans>0);
    figure, yyaxis left, histogram(panel_scan_Is,20), hold on;
    yyaxis right, plot(poi_x,poi_normy), hold off;
    title(['Panel#',num2str(panel),' for E range ',num2str(box_eV_range),' eV'])
end
%}    