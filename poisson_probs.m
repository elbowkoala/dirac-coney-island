 m%{
result1i = zeros(size(result1));
for i = 1:961
    ress = result1(:,:,i);
    ress = ress - 2.32;
    ress(ress<0) = 0;
    ress = round(ress ./ 4.21);
    result1i(:,:,i) = ress;
end
%}

BK = 5*K_bin*pix2invA;%.025 %+/- invA about k=0
BE1 = 36*E_bin*pix2eV; % bottom of box binding energy eV
BE2 = 30*E_bin*pix2eV; % top of box binding energy eV

disp({['BK=',num2str(BK)];[' BE1=',num2str(BE1)];[' BE2=',num2str(BE2)]})

bk = round(BK / pix2invA);
%disp(['bk/kbin=',num2str(bk/K_bin)])
be1 = round(BE1 / pix2eV);
be2 = round(BE2 / pix2eV);
boxI = 0;

indiv_PBs = zeros(1,length(find(involved_scans>0)));
indiv_SBI = zeros(1,length(find(involved_scans>0)));
ZZ = round(961*rand);


for iii = 1:size(region_list,1)
    nnn = 0;
    for jjj = find(involved_scans_panel(:,:,iii)>0)
        res1i = result1i(:,:,jjj);
        indiv_SBI(iii, jjj) = sum(sum(res1i(round(kLOS(jjj))-bk:round(kLOS(jjj))+bk,...
                                                                round(rfc_FL_Es(jjj))-be1:round(rfc_FL_Es(jjj))-be2)));
        panel_SBI(iii) = panel_SBI(iii) + indiv_SBI(iii,jjj);
        nnn = nnn+1;
    end
    panel_SBI(iii) = panel_SBI(iii) / nnn;
    
        



for i = find(involved_scans>0)
    res1i = result1i(:,:,i);
    indiv_SBI(i) = round(sum(sum(res1i(round(kLOS(i))-bk:round(kLOS(i))+bk,...
                                                        round(rfc_FL_Es(i))-be1:round(rfc_FL_Es(i))-be2))));     
    boxI = boxI + indiv_SBI(i);%round(sum(sum(res1i(round(kLOS(i))-bk:round(kLOS(i))+bk,...
                                                        %round(rfc_FL_Es(i))-be1:round(rfc_FL_Es(i))-be2))));                                                    
    %indiv_SBI(i) = round(sum(sum(res1i(round(kLOS(i))-bk:round(kLOS(i))+bk,...
                                                        %round(rfc_FL_Es(i))-be1:round(rfc_FL_Es(i))-be2))));
    %indiv_PBs(i) = (MBI^indiv_SBI(i))*exp(-MBI) / factorial(indiv_SBI(i));
     if i == ZZ
        figure, imagesc(imgaussfilt(res1i,5)), axis xy, hold on;
        plot([rfc_FL_Es(i),rfc_FL_Es(i)],[1,size(res1i,1)],'r'), hold on;
        plot([1,size(res1i,2)],[kLOS(i),kLOS(i)],'r'), hold on;
        plot([rfc_FL_Es(i)-be1,rfc_FL_Es(i)-be1,rfc_FL_Es(i)-be2,rfc_FL_Es(i)-be2,rfc_FL_Es(i)-be1],...
            [kLOS(i)-bk,kLOS(i)+bk,kLOS(i)+bk,kLOS(i)-bk,kLOS(i)-bk], 'r'), hold off;
        title(['i=',num2str(i)])
    end
end
%%Average events in range E_binding = [.20,.28] eV, k = [-.05, .05] invA
%%(per one scan)
mean_boxI = boxI / length(find(involved_scans>0));
MBI = mean_boxI;%(60+50+47+47+46+54+58+62+69)/9;%

for ii = find(involved_scans>0)
    res1i = result1i(:,:,ii);
    %boxI = boxI + sum(sum(res1i(round(kLOS(i))-bk:round(kLOS(i))+bk,...
                                                        %round(rfc_FL_Es(i))-be1:round(rfc_FL_Es(i))-be2)));
    %indiv_SBI(ii) = round(sum(sum(res1i(round(kLOS(ii))-bk:round(kLOS(ii))+bk,...
                                                        %round(rfc_FL_Es(ii))-be1:round(rfc_FL_Es(ii))-be2))));
    indiv_PBs(ii) = (MBI^indiv_SBI(ii))*exp(-MBI) / factorial(indiv_SBI(ii));
end



spec_boxIs = zeros(1,size(region_list,1));
figure,
II_n = 1;
for II = 1:size(region_list,1)
    spec = out_spec{1,II};
    eax = out_eax{1,II};
    kax = out_kax{1,II};
    box2 = find(eax<=-BE2,1,'last');%box1 + round((BE1-BE2)/pix2eV/E_bin);%
    box1 = box2 - round(abs(BE1-BE2)/pix2eV/E_bin);%find(eax<=-BE1,1,'last');
    K0p = round(size(spec,1)/2);
    
    box3 = K0p - round(bk/K_bin);% find(kax>-BK,1,'first')
    box4= K0p + round(bk/K_bin);%find(kax<+BK,1,'last')
    
    disp(['box size=',num2str(size(spec(box3:box4,box1:box2)))])

    
    spec_boxIs(II) = sum(sum(spec(box3:box4,box1:box2))) / out_nnn(II);
    
    SBI = round(spec_boxIs(II));
    
    poi = ((MBI^SBI)*(exp(-MBI))) / factorial(SBI)
    
    subplot(1,size(region_list,1),II_n)
    imagesc(kax, flip(eax), rot90(spec)), axis xy, hold on;
    plot([kax(K0p),kax(K0p)],[eax(1),eax(end)],'r'), hold on;
    plot([kax(box3),kax(box3),kax(box4),kax(box4),kax(box3)],...
        [eax(box1),eax(box2),eax(box2),eax(box1),eax(box1)],'r'), hold off;
    
    title({['SBI=',num2str(SBI)];['P=',num2str(poi)]})
    ylim([-.7,0.15])
    II_n = II_n+1;
    %{
    figure, 
    subplot(1,2,1)
    imagesc(kax, flip(eax), rot90(spec)), axis xy, hold on;
    plot([kax(K0p),kax(K0p)],[eax(1),eax(end)],'r'), hold on;
    plot([kax(box3),kax(box3),kax(box4),kax(box4),kax(box3)],...
        [eax(box1),eax(box2),eax(box2),eax(box1),eax(box1)],'r'), hold off;
    subplot(1,2,2)
    imagesc(kax(box3:box4),flip(eax(box1:box2)),rot90(spec(box3:box4,box1:box2))), axis xy, title(['II=',num2str(II)])
    %}
end
suptitle(['MBI=',num2str(MBI)])

poi_x = [1:1:100];
poi_yfactorials = zeros(size(poi_x));
for ii = poi_x
    poi_yfactorials(ii) = factorial(ii);
end
poi_y = ((MBI.^poi_x)*exp(-MBI)) ./ poi_yfactorials;
%figure, plot(poi_x, poi_y)

figure, histogram(indiv_SBI(indiv_SBI>0),50), hold on;
yyaxis right
plot(poi_x,poi_y)