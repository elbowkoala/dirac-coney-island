function [G] = GuassLowpass(imag,D0)%高斯低通滤波器模糊图像

F = fftshift(fft2(imag,size(imag,1),size(imag,2)));
u = 1:size(F,1);
v = 1:size(F,2);
[V,U] = meshgrid(v,u);
D = sqrt(( U -(floor(size(F,1)/2)+1)).^2 + ( V - (floor(size(F,2)/2)+1)).^2);
H = zeros(size(imag));
H = exp(-D.^2./(2*D0.^2));
G = F.*H;%guassian lowpass
G = ifft2(ifftshift(G));
end

