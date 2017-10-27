
%load rfc_FL_scan_170927.mat;
%load rfc_ncorr_scan_171019w.mat;

global x_x y_y FL_param
global x__x y__y %bcb_find dos_start

i = 946;%round(961*rand);
rfc_ncorr_scan

figure;

ss = imgaussfilt(cones(:,:,i),10);   
ssKpm = 5;
ssK = rfc_ks_after(i);
sss = sum( ss(ssK-ssKpm : ssK+ssKpm,:));
smin = find(sss(:,400:550)==min(sss(:,400:550)));
marj = 20;

y_y = sss(smin+399-marj:end);
x_x = [1:length(y_y)];
FL_param = rfc_FL_Es(i) - (smin+399-marj);

a1 = 40;
a2 = sss(smin+399);
a3 = max(y_y)-a2;
a4 = 0;

a0 = [a1, a2, a3, a4];
afinal = round(fminsearch(@bcbdevsum,a0));

SS_tot = sum((y_y-mean(y_y)).^2);
Rsquared = 1 - bcbdevsum(afinal)/SS_tot;

afinal(1) = round(abs(afinal(1)));
afinal(2) = abs(afinal(2));
afinal(3) = abs(afinal(3));
afinal(4) = abs(afinal(4));

y_fit = zeros(size(x_x));
y_fit(1:afinal(1)) = afinal(2);
y_fit(afinal(1)+1:round(FL_param)) = afinal(2)+afinal(3);
y_fit(round(FL_param)+1:end) = afinal(4);



kwts = abs([1:size(ss,1)]-ssK);
kwts(kwts>150) = 0;
DOS = sum(repmat(kwts, size(ss,2), 1)' .*ss);


bcb_find = smin + 399 - marj + afinal(1);
dos_start = find(DOS(150:250) == max(DOS(150:250))) + 149;
y__y = DOS(dos_start : bcb_find);
x__x = (1:length(y__y));

b1 = round(length(x__x)/2);
b2 = y__y(1)/b1;
b3 = y__y(1);
b4 = b2;
b5 = -b3;
b0 = [b1 b2 b3 b4 b5];
bfinal =  (fminsearch(@dosdevsum,b0));

bfinal(1) = round(abs(bfinal(1)));
bfinal(2) = -abs(bfinal(2));
bfinal(4) = +abs(bfinal(4));

y_fitb = zeros(size(x__x));
y_fitb(1:bfinal(1)) = bfinal(2)*x__x(1:bfinal(1)) + bfinal(3);
y_fitb(bfinal(1)+1:end) = bfinal(4)*x__x(bfinal(1)+1:end) + bfinal(5);

intercept = (bfinal(5) - bfinal(3))/(bfinal(2) - bfinal(4));
R__squared = 1 - (dosdevsum(bfinal))/SS__tot
step_E = intercept + dos_start - 1; 

subplot(2,2,1)
plot((1:length(sss)), sss, 'k'), hold on;%,(YY+i-400),'k'), hold on;
plot(smin+399,sss(smin+399),'r*'), hold on;
plot((x_x+smin+399-marj), y_fit, 'r'), hold off
xlim([1,size(ss,2)])
title(['Rsq=',num2str(Rsquared)])

subplot(2,2,2), plot(DOS)
title(['E=',num2str(step_E)])


subplot(2,2,3)
imagesc(imgaussfilt(cones(:,:,i),5)), axis xy, hold on;
plot([smin+399-marj+afinal(1),smin+399-marj+afinal(1)],[1,size(ss,1)],'w'), hold on;
plot(step_E,rfc_ks_after(i),'w*'), hold off
title(['DPI=',num2str(DPI_big(i))])

subplot(2,2,4)
plot(x__x, y__y,'k'), hold on;
plot(x__x, y_fitb, 'r'), hold off;
title(['Rsq=',num2str(R__squared)])

suptitle(['i=',num2str(i)])