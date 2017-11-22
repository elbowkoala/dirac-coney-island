%{
result1i = zeros(size(result1));
for i = 1:961
    ress = result1(:,:,i);
    ress = ress - 2.32;
    ress(ress<0) = 0;
    ress = round(ress ./ 4.21);
    result1i(:,:,i) = ress;
end
%}
BK = .023; %+/- invA about k=0
BE1 = 0.30; % bottom of box binding energy eV
BE2 = 0.25; % top of box binding energy eV

bk = round(BK / pix2invA);
disp(['bk/kbin=',num2str(bk/K_bin)])
be1 = round(BE1 / pix2eV);
be2 = round(BE2 / pix2eV);
boxI = 0;
for i = find(involved_scans>0)
    res1i = result1i(:,:,i);
    boxI = boxI + sum(sum(res1i(round(kLOS(i))-bk:round(kLOS(i))+bk,round(rfc_FL_Es(i))-be1:round(rfc_FL_Es(i))-be2)));
end
%%Average events in range E_binding = [.20,.28] eV, k = [-.05, .05] invA
%%(per one scan)
mean_boxI = boxI / length(find(involved_scans>0));
MBI = mean_boxI

spec_boxIs = zeros(1,size(region_list,1));
for II = 1:size(region_list,1)
    spec = out_spec{1,II};
    eax = out_eax{1,II};
    kax = out_kax{1,1};
    box1 = find(eax>-BE1,1,'first');
    box2 = find(eax>-BE2,1,'first');
    K0p = round(size(spec,1)/2);
    
    box3 = K0p - round(bk/K_bin);% find(kax>-BK,1,'first')
    box4= K0p + round(bk/K_bin);%find(kax<+BK,1,'last')
    
    disp(['box size=',num2str(size(spec(box3:box4,box1:box2)))])

    
    spec_boxIs(II) = sum(sum(spec(box3:box4,box1:box2))) / out_nnn(II);
    
    SBI = round(spec_boxIs(II));
    
    poi = ((MBI^SBI)*(exp(-MBI))) / factorial(SBI)
    
    figure, 
    subplot(1,2,1)
    imagesc(kax, flip(eax), rot90(spec)), axis xy
    subplot(1,2,2)
    imagesc(kax(box3:box4),flip(eax(box1:box2)),rot90(spec(box3:box4,box1:box2))), axis xy, title(['II=',num2str(II)])
end



