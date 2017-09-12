%% test contrast enhancement algorithm
%=========================
clear L2;clear L3;clear L4;clear L5;clear L6;clear L7;clear L8;
clear L02;clear L1;clear L;

L = im;
[N, M]=size(L);


m=15;n=15;
n1=fix(n/2);m1=fix(m/2);

a=L(1,1); b=L(1,M);
c=L(N,1); d=L(N,M);
%% we
for i=1:n1  
    for j=1:m1
        L1(i,j)=a;    
        L3(i,j)=b;    
        L6(i,j)=c;    
        L8(i,j)=d;
    end
end

L2 = L(1,1:M); 
L2 = repmat(L2, n1, 1);

L7 = L(N,1:M); 
L7 = repmat(L7, n1, 1);

L4 = L(1:N,1);
L4 = repmat(L4, 1, m1);

L5 = L(1:N,M);
L5 = repmat(L5, 1, m1);

L1=[L1;L4];  
L1=[L1;L6];  
L1=L1';

L2=[L2;L];  
L2=[L2;L7]; 
L2=L2';

L3=[L3;L5];  
L3=[L3;L8];  
L3=L3'; 

L1=[L1;L2];
L1=[L1;L3]; 
Lr=L1';

clear L2;clear L3;clear L4;clear L5;clear L6;clear L7;clear L8;
clear L02;clear L1;clear L;

%% qweqeq
%=========================
%====Определение среднеквадратических отклонений значений яркостей элементов локальной окрестности====
SV=zeros(N+2*n1,M+2*m1);
for i=1+n1:N+n1;
for j=1+m1:M+m1;
                 if j==1+m1;
                        D=0;
                        for a=-n1:n1;
                        for b=-m1:m1;
                           D(n1+1+a,m1+1+b)=Lr(i+a,j+b);
                        end;
                        end;
                 end;
           if j>1+m1;
            for a=-n1:n1;
              D(n1+1+a,m+1)=Lr(i+a,j+m1);
            end;
             D=D(1:n,2:m+1);
          end;            
               D=D(:);
               SV(i,j)=std(D);
 end;
 end;
%% qwe
%=====Определение массива усредненных значений элементов изображения=====
n_filter=3;m_filter=n_filter;
F=ones(n_filter,m_filter);
Lser=filter2(F,Lr,'same')/(n_filter*m_filter);
amin=.35;amax=.95;
C=(Lr-Lser)./(Lr+Lser+eps);
C=abs(C);
for i=1+n1:N+n1
 disp(i)
for j=1+m1:M+m1
      if SV(i,j)<=1;
        SV(i,j)=1;
      end;
                 if j==1+m1;
                        TM=0;
                        for a=-n1:n1;
                        for b=-m1:m1;
                           TM(n1+1+a,m1+1+b)=SV(i+a,j+b);
                        end;
                        end;
                 end;
         if j>1+m1;
            for a=-n1:n1;
              TM(n1+1+a,m+1)=SV(i+a,j+m1);
            end;
             TM=TM(1:n,2:m+1);
         end;    
    SV_MIN=min(min(TM));
    SV_MAX=max(max(TM));
            C(i,j)=C(i,j)^(amin+(amax-amin)*(SV(i,j)-SV_MIN)/(SV_MAX-SV_MIN));
    if Lroshyrena(i,j)>=Lser(i,j);
      Lvyh(i,j)=Lser(i,j)*(1+C(i,j))/(1-C(i,j));
    else
      Lvyh(i,j)=Lser(i,j)*(1-C(i,j))/(1+C(i,j));
    end;
   
%=====Проверка корректности диапазона=====
   if Lvyh(i,j)>=255;      Lvyh(i,j)=255;   end;
   if Lvyh(i,j)<=0;      Lvyh(i,j)=0;   end;
end;
end;
Lvyh=Lvyh(n1+1:N+n1,m1+1:M+m1);
Lvyh=round(Lvyh);
L=Lroshyrena(n1+1:N+n1,m1+1:M+m1);

%% wrw
%=====Визуализация результатов=====
colormap(gray(255));
subplot(221);image(L);axis('image');
subplot(222);image(Lvyh);axis('image');
