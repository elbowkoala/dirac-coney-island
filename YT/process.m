clear all
close all 
clc

load cones.mat
i=1,j=1;
m=1,n=1,k=1;
B=zeros(300,800,961);
C=zeros(1,10000000);
global A

result=zeros(300,800,961);


for k=1:961
if rem(k,50) == 0
    disp(['on YT scan ',num2str(k)])
end   
    A=zeros(302,802);
    
    for m=1:300
        for n=1:800
            if cones(m,n,k)~=0
                A(m+1,n+1)=1;
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
        result(m,n,k)=A(m+1,n+1);
    end
end


end

test=zeros(300,800);
for m=1:300
    for n=1:800
        test(m,n)=result(m,n,961);
    end
end

mesh(test)
%mesh(ek)
