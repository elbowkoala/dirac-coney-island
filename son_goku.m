goku = zeros(size(draw_scan_window));
%figure, imagesc(goku), axis xy
draw_box = zeros(size(fass));

LWP = 2;
LWN = 3;

A = 1.9;
B = 0; 

E_0 = 70;

K0 = round(size(draw_box,1)/2);
K_0s      = [K0-1, K0, K0+2, K0+1, K0+7, K0+8, K0+10, K0+11, K0+13];
I_factors = [2,  1, 2,    1.5,    .5,    2,    2,     1,     1, 1];
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
krillins = imnoise(krillins,'salt & pepper',.2);
figure, imagesc(imgaussfilt(krillins,1)), axis xy, hold on;
plot([E_0,E_0],[1,size(draw_box,1)],'r'), hold off

    
 %{   
K_1 = K_0 + 7;
draw_x0 = [(1:size(draw_box,1))-K_0]';
draw_x1 = [(1:size(draw_box,1))-K_1]';
draw_yp = A*(draw_x0) + B*((draw_x0).^2) + E_0;
draw_yn = -A*(draw_x0) + B*((draw_x0).^2) + E_0;
draw_yp1 = A*(draw_x1) + B*((draw_x1).^2) + E_0;
draw_yn1 = -A*(draw_x1) + B*((draw_x1).^2) + E_0;

draw_itp = horzcat(draw_yp,draw_x0+K_0);
draw_itn = horzcat(draw_yn,draw_x0+K_0);
draw_itp1 = horzcat(draw_yp1,draw_x1+K_1);
draw_itn1 = horzcat(draw_yn1,draw_x1+K_1);

draw_itp_curt = draw_itp(max(1,round(-A/(2*B))+K_0):end,:);
draw_itn_curt = draw_itn(1:min(length(draw_x0),round(K_0+A/(2*B))),:);

draw_itp_it = reshape(draw_itp_curt',1,[]);
draw_itn_it = reshape(draw_itn_curt',1,[]);
draw_itn1_it = reshape(draw_itn1',1,[]);
draw_itp1_it = reshape(draw_itp1',1,[]);
%if length(B_range) == 1
    draw_itp_it_ = zeros(1,4);
    draw_itn_it_ = zeros(1,4);
    draw_itp1_it_ = zeros(1,4);
    draw_itn1_it_ = zeros(1,4);
    draw_itp_it_(1:2) = draw_itp_it(1:2);
    draw_itn_it_(1:2) = draw_itn_it(1:2);
    draw_itp1_it_(1:2) = draw_itp1_it(1:2);
    draw_itn1_it_(1:2) = draw_itn1_it(1:2);
    draw_itp_it_(3:4) = draw_itp_it(end-1:end);
    draw_itn_it_(3:4) = draw_itn_it(end-1:end);
    draw_itp1_it_(3:4) = draw_itp1_it(end-1:end);
    draw_itn1_it_(3:4) = draw_itn1_it(end-1:end);
    draw_itp_it = draw_itp_it_;
    draw_itn_it = draw_itn_it_;
    draw_itp1_it = draw_itp1_it_;
    draw_itn1_it = draw_itn1_it_;
%end      
ITP = insertShape(draw_box, 'Line', draw_itp_it,'LineWidth',LWP);
ITP1 = insertShape(draw_box,'Line', draw_itp1_it, 'LineWidth',LWP1);
%{
if nnz(ITP(:,1))==0   %Only consider A,B values that give lines that extend across draw_box
    %ABBA_norm_table(A_i,B_i) = NaN;
    ABBA_IT_i = ABBA_IT_i + 1;
    continue
end
%}
ITN = insertShape(draw_box, 'Line', draw_itn_it, 'LineWidth', LWN);
ITN1 = insertShape(draw_box, 'Line', draw_itn1_it, 'LineWidth',LWN1);
ITP = ITP(:,:,1);
ITP1 = ITP1(:,:,1);
ITN = ITN(:,:,1);
ITN1 = ITN1(:,:,1);

krillin = ITP+ITN;
krillin1 = ITP1+ITN1;

figure, imagesc(krillin+krillin1), axis xy, hold on;
plot([E_0,E_0],[1,size(draw_box,1)],'r'), hold off

krillin(krillin~=0)=1;
krillin = abs(1-krillin);

gohan = bwconncomp(krillin);
krillin(gohan.PixelIdxList{1})=0;
krillin(gohan.PixelIdxList{2})=0;
krillin(gohan.PixelIdxList{3})=0;

%figure, imagesc(krillin), axis xy




%figure, 
%imagesc(imgaussfilt(ITP+ITN,2)), axis xy;
%}