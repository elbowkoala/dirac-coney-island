clear all
close all 
clc

load H.mat
cutoffintensity=35;
cutoffbins=round(cutoffintensity/0.064466);

l=linspace(1,10000,10000);
Z=l*H';

sum=0;
for i=1:cutoffbins
    sum=s/um+l(1,i)*H(1,i);
end

sum/Z
