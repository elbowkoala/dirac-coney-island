
i = 34;
cone = cones(:,:,i);
coneb = Binning_2d(cone,bin_E,bin_k);
conebf = imgaussfilt(coneb,fass_sigma); 
corn = conebf;

figure;
for MDC = 10 : 5 : size(corn,2)-10 
    strip = corn(:,MDC-1:MDC+1);
    Strip = sum(strip,2);
    plot([1:length(Strip)], Strip + MDC), hold on;
end
