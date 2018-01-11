clear all
close all 
clc

load cones.mat
i=1,j=1;
m=1,n=1,k=1;
B=zeros(300,800,961);
C=zeros(1,10000000);
global A
global A1
global A2
global I0
global p

I0=zeros(1,1000000);
p=0;

result1=zeros(300,800,961);
result2=zeros(300,800,961);

for k=1:961
    
A=zeros(302,802);
A1=zeros(302,802);
A2=zeros(302,802);

for i=2:301
 for j=2:801
     a=rand;
     if a>0.5
  A(i,j)=a;
     end
 end
end

for m=1:300
    for n=1:800
        if A(m+1,n+1)~=0
            Search1(m+1,n+1);
        end
    end
end

for m=1:300
    for n=1:800
        result1(m,n,k)=A1(m+1,n+1);
        result2(m,n,k)=A2(m+1,n+1);
    end
end

end
mesh(A1)




