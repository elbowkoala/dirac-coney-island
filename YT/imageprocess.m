clear all
close all 
clc

load cones5.mat

Z=zeros(300,800);

for i=1:961
Z=Z+result1(:,:,i);
end

mesh(Z)