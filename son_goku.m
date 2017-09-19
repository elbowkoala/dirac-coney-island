%goku = zeros(size(draw_scan_window));
%figure, imagesc(goku), axis xy
%draw_box = zeros(size(fass));

LWP = 1;
LWN = 1;

A = 1.9;
B = 0; 



parab_A = 0.12;
parab_box = zeros(size(fass,1),90);
parab_x = [1:size(parab_box,1)]';
parab_K0 = round(length(parab_x)/2);

parab_E_off = 50;
cond_band_window = fass(:,parab_E_off+(1:size(parab_box,2)));
parab_E0 = 60;%size(parab_box,2)-size(cond_band_window,2);
parab = parab_A * (parab_x-parab_K0).^2 + parab_E0;

parabb = horzcat(parab,parab_x);
parabb = reshape(parabb',1,[]);
cond_band = insertShape(parab_box,'Line',parabb);
cond_band = cond_band(:,:,1);
figure, 
%subplot(1,2,1), imagesc(cond_band), axis xy,
%subplot(1,2,2), imagesc(cond_band_window), axis xy

E = mat2gray(norman(cond_band_window,0,3));
Ea = rot90(E,-1);
I = mat2gray(cond_band);
Ia = rot90(I,-1);
greenlol = cat(3,ones(size(Ea)),zeros(size(Ea)),zeros(size(Ea)));
imshow(Ea,'InitialMag','fit','DisplayRange',[0 1]), hold on;
h = imshow(greenlol);  axis xy
set(h,'AlphaData',Ia);
hold off;

%{
K0 = round(size(draw_box,1)/2);
K_0s     = [K0, K0+20,K0+21,K0+22];
I_factors = [2,  1, 1, 1];%, 2,    1.5,    .5,    2,    2,     1,     1, 1];

krillins = zeros(size(draw_box));
for K_0_i = 1:length(K_0s)
    K_0 = K_0s(K_0_i);
    
    draw_x0 = [(1:size(draw_box,1))-K_0]';
    draw_yp = A*(draw_x0) + B*((draw_x0).^2) + E_0;
    draw_yn = -A*(draw_x0) + B*((draw_x0).^2) + E_0;
    draw_itp = horzcat(draw_yp,draw_x0+K_0);
    draw_itn = horzcat(draw_yn,draw_x0+K_0);
    draw_itp_curt = draw_itp(max(1,round(-A/(2*B))+K_0):end,:);
    draw_itn_curt = draw_itn(1:min(length(draw_x0),round(K_0+A/(2*B))),:);

    draw_itp_it = reshape(draw_itp_curt',1,[]);
    draw_itn_it = reshape(draw_itn_curt',1,[]);
    %if length(B_range) == 1
    draw_itp_it_ = zeros(1,4);
    draw_itn_it_ = zeros(1,4);
    draw_itp_it_(1:2) = draw_itp_it(1:2);
    draw_itn_it_(1:2) = draw_itn_it(1:2);
    draw_itp_it_(3:4) = draw_itp_it(end-1:end);
    draw_itn_it_(3:4) = draw_itn_it(end-1:end);
    draw_itp_it = draw_itp_it_;
    draw_itn_it = draw_itn_it_;
    %end      
    ITP = insertShape(draw_box, 'Line', draw_itp_it,'LineWidth',LWP);
    ITN = insertShape(draw_box, 'Line', draw_itn_it, 'LineWidth', LWN);
    ITP = ITP(:,:,1);
    ITN = ITN(:,:,1);
    krillin = ITP+ITN;
    krillin = krillin.*I_factors(K_0_i);
    krillins = krillins+krillin;
end

%krillins = imnoise(krillins,'salt & pepper',.2);

figure, imagesc(imgaussfilt(krillins,1)), axis xy, hold on;
plot([E_0,E_0],[1,size(draw_box,1)],'r'), hold off
%}


