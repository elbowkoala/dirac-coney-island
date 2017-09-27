

%fermi_= zeros(1,num_scans);


for i=1:num_scans

    duhh = sum(imgaussfilt(cones(:,:,i),6),1);
    duhhh = smooth(duhh,3,'rlowess');
    duhhh = duhhh./max(duhhh);
    duh_bkgd =  mean(duhhh(700:780));
    if rem(i,100)==0;
        disp(i)
    end
    fermi_(i) = find(duhhh(550:700)>=0.5+duh_bkgd/2,1,'last') + 549;
    %{
    figure, subplot(2,1,1), plot(duhhh), hold on;title(num2str(i))
    plot(fermi_(i),duhhh(fermi_(i)),'r*'), hold off
    subplot(2,1,2), imagesc(imgaussfilt(cones(:,:,i),5)), axis xy, hold on;
    plot([fermi_(i),fermi_(i)],[1,300],'r'), hold off;
    title(num2str(fermi_(i)))
    pause(.001)
    %}
end

figure, 
plot(fermi_)
%subplot(3,1,1), plot(duhmaxAs,'r');
%subplot(3,1,2), plot(duhmaxBs,'b');
%subplot(3,1,3), plot(duhminCs,'p');