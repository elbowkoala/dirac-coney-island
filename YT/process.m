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
    
    for m=1:300
        for n=1:800
                A(m+1,n+1)=cones(m,n,k);          
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

test1=zeros(300,800);
test2=zeros(300,800);
for m=1:300
    for n=1:800
        test1(m,n)=result1(m,n,961);
        test2(m,n)=result2(m,n,961);
    end
end
I1=zeros(1,p);
for i=1:p
   I1(1,i)=I0(1,i);
end

% mesh(test1)

H=hist(I1,10000)

%mesh(test2)
%mesh(ek)
