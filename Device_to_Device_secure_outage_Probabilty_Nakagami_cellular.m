%
% Copyright (c) 2022, Ajay Kumar, and  Sudhan Majhi, IISc Bangalore.
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice, this
%   list of conditions and the following disclaimer.
% 2. Redistributions in binary form must reproduce the above copyright notice,
%   this list of conditions and the following disclaimer in the documentation
%   and/or other materials provided with the distribution.
% 3. The reference listed below should be cited if the corresponding codes are used for
%   publication..
%
%THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
%ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
%WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
%DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
%ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

%    - Freely distributed for educational and research purposes
%References

%  [R1]. A. Kumar, S. Majhi and H. -C. Wu, "Physical-Layer Security of Underlay MIMO-D2D Communications by Null Steering Method Over Nakagami-m and Norton Fading Channels," in IEEE Transactions on Wireless Communications, vol. 21, no. 11, pp. 9700-9711, Nov. 2022, doi: 10.1109/TWC.2022.3178758.
%  [R2]. A. Kumar, S. Perveen, S. Singh, A. Kumar, S. Majhi and S. K. Das, "6th Generation: Communication, Signal Processing, Advanced Infrastructure, Emerging Technologies and Challenges," 2021 6th International Conference on Computing, Communication and Security (ICCCS), Las Vegas, NV, USA, 2021, pp. 1-16, doi: 10.1109/ICCCS51487.2021.9776334.



clc;
clear all;
close all;

mc=1000; % monte carlo iterations
N=2; % Number of antenna at each Transmitter

% channel parameters (Nakagami)
m=[1 2 3 4]; % shape parameter
w=[1 2 1 3]; % scale parameter

R=2.5; % target data rate

%% channel intitialization
H11=[]; 
H12=[];
H13=[];
H1cu=[];
H21=[];
H22=[];
H23=[];
H2cu=[];
H31=[];
H32=[];
H33=[];
H3cu=[];
HB1=[];
HB2=[];
HB3=[];
HBcu=[];
H1a=[];
H21=[];
H3a=[];
Hcua=[];
%%
SNR_dB=0:1:25; % SNR in dB
outcu_wobf=zeros(1,length(SNR_dB)); % intialization of outage probability for CU
outcu=zeros(1,length(SNR_dB)); % intialization of outage probability for CU

n=1;
for x=SNR_dB
    x
    P0=10.^(x/10); %total transmit power in watts
    pcu=P0/2; %power for CU
    pd=P0/6;% power for DU
    c4_cu=0;%total secrecy for DU
    c4_cu_wobf=0;%total secrecy for DU
  
    for i=1:mc
        %% Nakagami channel generation
        for j=1:N    
        H11(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        H12(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        H13(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        H1cu(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        H21(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        H21(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        H23(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        H2cu(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        H31(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        H32(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        H33(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        H3cu(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        HB1(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        HB2(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        HB3(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        HBcu(j)=sqrt(gamrnd(m(j),w(j),[1,1]));
        end
        %%
        H1a=[H12;H13;H1cu]; 
        H2a=[H21;H23;H2cu];
        H3a=[H31;H32;H3cu];
        Hcua=[HB1;HB2;HB3];

        %% Null steering
        w1=null(H1a); % ZFBF or null steering
        w2=null(H2a); % ZFBF or null steering
        w3=null(H3a); % ZFBF or null steering
        wcu=null(Hcua); % ZFBF or null steering
       %% secrecy capacity and secrecy outage probability calucations
        S_cu_wobf=(pcu*sum(HBcu* HBcu'))/(pd*sum(H1cu*H1cu')+pd*sum(H2cu*H2cu')+pd*sum(H3cu*H3cu')+1); % received SNR at CU
        c3_cu_wobf=log2(1+S_cu_wobf); % secrecy capacity of CU
        c4_cu_wobf=c3_cu_wobf+c4_cu_wobf;

        S_cu=(pcu*sum( HBcu*HBcu'))/(pd*sum(H1cu*w1*w1'*H2cu')+pd*sum(H2cu*w2*w2'*H2cu')+pd*sum(H3cu*w3*w3'*H3cu')+1); % received SNR at CU
        c3_cu=log2(1+S_cu); % secrecy capacity of CU
        c4_cu=c3_cu+c4_cu;
        % secrecy outage probability calculations
        if c3_cu_wobf<R
             outcu_wobf(n)=outcu_wobf(n)+1;
         end
         if c3_cu<R
             outcu(n)=outcu(n)+1;
         end
    end
     c4_cu11_wobf(n)=c4_cu_wobf/mc; 
     c4_cu11(n)=c4_cu/mc; 

        n=n+1;
end
plot(SNR_dB,c4_cu11_wobf,'LineWidth',2);
hold on;
plot(SNR_dB,c4_cu11,'LineWidth',2);
xlabel('SNR (dB)');
ylabel('Total secrecy capacity');
legend('without null steering','with null steering')
title('Total secrecy capacity');
figure;
semilogy(SNR_dB,outcu_wobf/mc,'LineWidth',2);
hold on;
semilogy(SNR_dB,outcu/mc,'LineWidth',2);
xlabel('SNR (dB)');
ylabel('CU secrecy outage probability');
legend('without null steering','with null steering')
title('Total secrecy outage probability for CU');

