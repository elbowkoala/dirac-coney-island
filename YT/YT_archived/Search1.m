function[]=Search1(i,j)
global A
k=zeros(2,500);
k(1,1)=i;%��ǰ�Ľڵ�the present position of searching(k(1,1),k(2,1))
k(2,1)=j;%��ǰ�Ľڵ�the present position of searching(k(1,1),k(2,1))
n=1;
conti=1;%1 means continue,0 means stop

A(i,j)=0;

i=0;
j=0;%�����������ѭ����ʼ������ע�������Ͻǵ����ز����Ƿ���ģ�����������
    %This is to set the loop start, but beware that the upper left square
    %of the whole picture couldn't be non zero

while conti==1
%     i~=k(1,1) | j~=k(2,1)

   i=k(1,n);j=k(2,n);%�����Ϊ���ڵ�һ�ν���ѭ����ʱ��ij��ȡֵ����
            %to set the value of i,j to be normal
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
 if A(i-1,j)~=0
       A(i-1,j)=0;
       i=i-1;
       n=n+1;
       k(1,n)=i;
       k(2,n)=j;
       conti=1;
       
   else if A(i,j+1)~=0
           A(i,j+1)=0;
           j=j+1;
           n=n+1;
           k(1,n)=i;
           k(2,n)=j;
           conti=1;
           
       else if A(i+1,j)~=0
               A(i+1,j)=0;
               i=i+1;
               n=n+1;
               k(1,n)=i;
               k(2,n)=j;
               conti=1;
               
           else if A(i,j-1)~=0
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

%%%%%%%%The while���������ends here%%%%%%%%%%%%%%%%%%
   
   A(k(1,1),k(2,1))=1;
       
    
end