global XX
global YY

%asdf = imgaussfilt(cones(:,:,round(961*rand)),7);
asdf = imgaussfilt(regional_arpes(:,:,3),5)+imgaussfilt(regional_arpes(:,:,4),5);
%asdf = imgaussfilt(raw_full_cone,5);
interped = 1;

if interped == 1
    MDC_steps = [400:10:630];
    XX = (K_interp(:,1).*pix2invA)';
    a2 = -.05; a5 = +.05;
    shift_factor = 1;
end

if interped == 0
    MDC_steps = [400:5:500];
    XX = (1:size(asdf,1)).*pix2invA';
    a2 = XX(250);  a5 = XX(100);
    shift_factor = 5;
end

figure;
subplot(1,2,1),
lorentz_params = cell(length(MDC_steps),6);
lorentz_params(1,:) = {'a1/(2pi)','a2','2*sqrt(a3)','a4/(2pi)','a5','2*sqrt(a6)'};
step_n = 2;
for i = MDC_steps
  
    asdff = asdf(:,i);    
    YY = asdff';
    
    a3 = ((max(XX) - min(XX))/10)^2;
    a1 = max(YY)*a3;

    a6 = a3;
    a4 = a1;

    a0 = [a1, a2, a3, a4, a5, a6];
    afinal = fminsearch(@trydevsum,a0);
    afinal(1) = abs(afinal(1));
    afinal(4) = abs(afinal(4));
    afinal(3) = abs(afinal(3));
    afinal(6) = abs(afinal(6));
    
    lorentz_params(step_n,:) = {afinal(1)/(2*pi),afinal(2), 2*sqrt(afinal(3)), afinal(4)/(2*pi),afinal(5), sqrt(afinal(6))};
    step_n = step_n + 1;
    
    %plot(XX,(YY+i-400),'k'), hold on;
    plot(XX, ((i-MDC_steps(1))*shift_factor+abs(afinal(1))./((XX-afinal(2)).^2+afinal(3)) + abs(afinal(4))./((XX-afinal(5)).^2+afinal(6))),'r'), hold on;
end
subplot(1,2,2), imagesc(asdf'), axis xy
lorentz_params