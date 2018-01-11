clear all
close all 
clc

load finalresult.mat
 h=histogram(I1,10000);
 
 Hx=linspace(0,644.66,10000);
 
 H=[h.Values];
 
 saddles=round([4.126 8.187 12.51 16.95 21.47 25.66 29.72 33.97 38.1 42.48 46.61 50.86 55.12 59.18 63.31 67.56]/0.064466);
 
 for i=1:15
     weight=0;
     norm=0;
     for j=saddles(1,i):saddles(1,i+1)
         weight=weight+H(1,j)*j;
         norm=norm+H(1,j);
     end
     COMs(1,i)=round(weight/norm);
 end
 
 COMs=COMs*0.064466
 
 histogram(I1,10000);
 hold on
 for i=1:15
 plot([COMs(1,i),COMs(1,i)],[0,10000],'R');
 hold on
 end
 