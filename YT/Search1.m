function[]=Search1(i,j)
global A
global A1
global A2
global I0
global p

p=p+1;

k=zeros(2,500);
k(1,1)=i;%当前的节点the present position of searching(k(1,1),k(2,1))
k(2,1)=j;%当前的节点the present position of searching(k(1,1),k(2,1))
n=1;
conti=1;%1 means continue,0 means stop
I=0;%亮度为零 Intensity is set to be zero in the beginning
iweigh=0;%i坐标的权重 the weigh of i
jweigh=0;%
ibar=0;%i的平均坐标 center of mass of x axis
jbar=0;%


%A(i,j)=0;

i=0;
j=0;%这个设置是让循环开始，所以注意最左上角的像素不能是非零的，否则会出问题
    %This is to set the loop start, but beware that the upper left square
    %of the whole picture couldn't be non zero

while conti==1
%     i~=k(1,1) | j~=k(2,1)

   i=k(1,n);j=k(2,n);%这个是为了在第一次进入循环的时候ij的取值正常
            %to set the value of i,j to be normal
   I=I+A(i,j);
   iweigh=iweigh+i*A(i,j);
   jweigh=jweigh+j*A(i,j);
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
 if A(i-1,j)~=0
       I=I+A(i-1,j);
       iweigh=iweigh+i*A(i-1,j);
       jweigh=jweigh+j*A(i-1,j);
       A(i-1,j)=0;
       i=i-1;
       n=n+1;
       k(1,n)=i;
       k(2,n)=j;
       conti=1;
       
   else if A(i,j+1)~=0
           I=I+A(i,j+1);
   iweigh=iweigh+i*A(i,j+1);
   jweigh=jweigh+j*A(i,j+1);
           A(i,j+1)=0;
           j=j+1;
           n=n+1;
           k(1,n)=i;
           k(2,n)=j;
           conti=1;
           
       else if A(i+1,j)~=0
               I=I+A(i+1,j);
   iweigh=iweigh+i*A(i+1,j);
   jweigh=jweigh+j*A(i+1,j);
               A(i+1,j)=0;
               i=i+1;
               n=n+1;
               k(1,n)=i;
               k(2,n)=j;
               conti=1;
               
           else if A(i,j-1)~=0
                   I=I+A(i,j-1);
   iweigh=iweigh+i*A(i,j-1);
   jweigh=jweigh+j*A(i,j-1);
                   A(i,j-1)=0;
                   j=j-1;
                   n=n+1;
                   k(1,n)=i;
                   k(2,n)=j;
                   conti=1;
                   
               else
                   n=n-1;
                   if n==0
                       conti=0;
                   else
                       i=k(1,n);j=k(2,n);
                       conti=1;
                   end
                  
                   
                   end  
               end
           end
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
end

%%%%%%%%The while到这里结束ends here%%%%%%%%%%%%%%%%%%
   
   I;i;j;iweigh;jweigh;
% if I>0.2

% saddles=[4.126 8.187 12.51 16.95 21.47 25.66 29.72 33.97 38.1 42.48 46.61 50.86 55.12 59.18 63.31 67.56];
% COMs =[6.3177 10.6369 14.9561 19.2109 23.5301 27.7204 31.8462 36.0365 40.2268 44.4815 48.6718 52.8621 57.0524 61.1782 65.3685];

   ibar=round(iweigh/I);
   jbar=round(jweigh/I);
% for i=1:15
%     if I>saddles(1,i) && I<saddles(1,i+1)
%    I=COMs(1,i);
%     end
% end
%    I;
   A1(ibar,jbar)=I;
   A2(ibar,jbar)=1;
   
   % end    
   %A(k(1,1),k(2,1))=1;
   
   I0(1,p)=I;
   
end