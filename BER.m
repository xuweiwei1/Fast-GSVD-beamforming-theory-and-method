clc;
clear;
close all;
warning off;
addpath(genpath(pwd));
rng('default')
%FPGA/MATLAB/simulink
Nt =4;
Nri1 =2;
Nri2 = 2;
S=4;
Tc = 1100;
sigma2=1;
SNRdbs=-10:2:30;
BER1=zeros(S,length(SNRdbs));
BER2=zeros(S,length(SNRdbs));
t1=0;

for i_SNR=1:length(SNRdbs)
count1=zeros(1,S);
count2=zeros(1,S);
        SNR=10^(SNRdbs(1,i_SNR)/10);
        P=sigma2*SNR;
        loop_num=1000;      
 for i_loop=1:loop_num
            q=rand(S,Tc)>=0.5;
            s=QPSK_mapper(q,Tc,S,1);
            H1_real=randn(Nri1,Nt);
            H1_imag=randn(Nri1,Nt);
            H1=(1/sqrt(2))*complex(H1_real,H1_imag);
            H2_real=randn(Nri2,Nt);
            H2_imag=randn(Nri2,Nt);
            H2=(1/sqrt(2))*complex(H2_real,H2_imag);
            n1_real=sqrt(sigma2/2)*randn(Nri1,Tc/2);
            n1_imag=sqrt(sigma2/2)*randn(Nri1,Tc/2);
            n1=complex(n1_real,n1_imag);
            n2_real=sqrt(sigma2/2)*randn(Nri2,Tc/2);
            n2_imag=sqrt(sigma2/2)*randn(Nri2,Tc/2);
            n2=complex(n2_real,n2_imag);
            [U,V,X,~,~]=fgsvd(H1,H2);
            [u,v,c]=svd(X);
            pre=c*pinv(v)*u';  
            y1=sqrt(P/Nt)*H1*pre*s+n1;
            yy1=U'*y1;
            y2=sqrt(P/Nt)*H2*pre*s+n2;
            yy2=V'*y2;        
            yy1=pinv(U'*sqrt(P/Nt)*H1*pre)*yy1;
            yy2=pinv(V'*sqrt(P/Nt)*H2*pre)*yy2;    
            qq1=zeros(S,Tc);
for i1=1:S
    for i2=1:Tc/2
        if real(yy1(i1,i2))>=0
            qq1(i1,2*i2-1)=0;
        else
            qq1(i1,2*i2-1)=1;
        end
        if imag(yy1(i1,i2))>=0
            qq1(i1,2*i2)=0;
        else
            qq1(i1,2*i2)=1;
        end
    end
end
qq2=zeros(S,Tc);
for i1=1:S
    for i2=1:Tc/2
        if real(yy2(i1,i2))>=0
            qq2(i1,2*i2-1)=0;
        else
            qq2(i1,2*i2-1)=1;
        end
        if imag(yy2(i1,i2))>=0
            qq2(i1,2*i2)=0;
        else
            qq2(i1,2*i2)=1;
        end
    end
end
for i1=1:S
    for i2=1:Tc
        if q(i1,i2)~=qq1(i1,i2)
            count1(1,i1)=count1(1,i1)+1;
        end
    end
end
for i1=1:S
    for i2=1:Tc
        if q(i1,i2)~=qq2(i1,i2)
            count2(1,i1)=count2(1,i1)+1;
        end
    end         
end
for i1=1:S
    BER1(i1,i_SNR)=count1(1,i1)/(loop_num*Tc);
    BER2(i1,i_SNR)=count2(1,i1)/(loop_num*Tc);
end
end
end
    figure(1)
    semilogy(SNRdbs,BER1(1,:),'b--','linewidth',1);
    hold on;
    semilogy(SNRdbs,BER1(2,:),'m:','linewidth',1);
    hold on;
   scatter(SNRdbs,BER2(3,:),'ro','filled');
    grid on;
    scatter(SNRdbs,BER2(4,:),'bd','filled');
    grid on;
legend('U_1,PC_1','U_1,PC_2','U_2,PC_1','U_2,PC_2');    
xlabel('SNR (dB)');
ylabel('symbol error rate');


