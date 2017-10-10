draw_box = zeros([71+10], [75+10]);
draw_x = (1:size(draw_box,1))';
draw_x0 = draw_x - round(size(draw_box,1)/2);

E_0 = 30;
K_0 = round(size(draw_box,1)/2);

A = 1.5;
B = .02;

draw_yp = A*(draw_x0) + B*((draw_x0).^2) + E_0;
draw_yn = -A*(draw_x0) + B*((draw_x0).^2) + E_0;
draw_itp = horzcat(draw_yp,draw_x);
draw_itn = horzcat(draw_yn,draw_x);

draw_itp_it = reshape(draw_itp',1,[]);
draw_itn_it = reshape(draw_itn',1,[]);

ITP = insertShape(draw_box,'Line',draw_itp_it);
ITN = insertShape(draw_box,'Line',draw_itn_it);

ITT = ITP(:,:,1)+ITN(:,:,1);
IT = ITT(6:end-5,6:end-5);
figure, 
imagesc(ITT), axis xy
figure, imagesc(IT), axis xy
