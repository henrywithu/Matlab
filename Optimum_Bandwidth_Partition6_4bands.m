clear;clc
%% Channel Model
Fs=800e6;  
Num=4096;
f=[1:1:400]*1e6;
l=0.025; %% unit is km
%l=0.150;
R0c=176.6*l;
Ac=0.0500079494;
L0=1090.8e-6*l;
L_inf=504.5e-6*l;
B=0.705;
Fm=32570;
C_inf=48.55e-9*l;
G0=1.47653e-9*l;
Ge=0.91;
R_f=1./(1./(R0c^4+Ac*f.^2).^0.25);
L_f=(L0+L_inf*(f/Fm).^B)./(1+(f/Fm).^B);
C_f=C_inf;
G_f=G0*f.^Ge;
Gama=((R_f+i*2*pi*f.*L_f).*(G_f+i*2*pi*f.*C_f)).^0.5;
Z0=((R_f+i*2*pi*f.*L_f)./(G_f+i*2*pi*f.*C_f)).^0.5;
Zl=Z0;                                              %% perfect load impedance matching
Z1=(Zl+Z0.*tanh(Gama*l))./(1+Zl./Z0.*tanh(Gama*l));  %% input impedance of channel with Zl load
Zs=Z1;                                              %% perfect source impedance matching
H_f=(Z1+Zs)./Z1.*Z0.*sech(Gama*l)./((Zs.*(Z0./Zl+tanh(Gama*l)))+(Z0.*(1+Z0./Zl.*tanh(Gama*l))));   %% transfer function of channel
H_f_amp=20*log10(abs(H_f));
%figure(1);
[b,a] = invfreqz(abs(H_f),2*pi*f./Fs,20,20);    %% transfer H_f to H_Z
%[b,a] = invfreqz(H_f,2*pi*f./Fs,30,30,[],300000);
[h,w] = freqz(b,a,Num,Fs);
%[h2,t] = impz(b,a);
% p_sq=sum(h2.^2)
% p_sq2=sum(b.^2)/sum(a.^2)
 figure(1);
 plot(f/1e6,H_f_amp,'r-','LineWidth',2);
 hold on; 
 plot(w/1e6,20*log10(abs(h)),'k*','LineWidth',2);
 grid on;
 xlabel('f(MHz) ','fontsize',14);
 ylabel('Attenuation(dB)','fontsize',14);
% figure(2);
% plot(t,h2,'rd','LineWidth',2);
% break;
%% Bandwidth partition
b0=15;
bs=1;
Band_total=400;  %% MHz
Fs2=Band_total*2;      %% MHz, for FREQZ transform 
PARr=10^(15.2/10);  %% Receiver peak to average ratio
Epsilon_x=10^(-60/10);     %% transimited signal energy
Sigma_r=10^(-130/10);  %% total receiver noise
Gap=10^(7.2/10);  %% SNR Gap

ai=0.0;  %$ excess bandwidth 
gi=0.2;    %% guard band

Nc=32;           %% Number of DMT sub carrier
Nband=2;         %% Number of partitioned channel 

fs_init=[50;50;50;50];
p_sq2=sum(b.^2)/sum(a.^2)
%p_sq=sum(h2.^2)
[H_freq,f2]= freqz(b,a,Num,Fs2*1e6);


%% Caculation for fs_init

sub_band=2*Num*fs_init/Fs2/Nc;                                                 %% bandwidth of sub carrier in DMT
%offset_2=[ones(1,Nc)*Num/Fs2*fs_init(1)*gi; ones(1,Nc)*(2*Num/Fs2*fs_init(1)*(1+gi)+Num/Fs2*fs_init(2)*gi)];   %% half guard band before the channel and other half guard after the channel. don't consider the excess bandwidth  

offset_4=[ones(1,Nc)*Num/Fs2*fs_init(1)*0.5*gi;  ones(1,Nc)*(2*Num/Fs2*fs_init(1)*(1+gi)+Num/Fs2*fs_init(2)*0.5*gi);  ones(1,Nc)*(2*Num/Fs2*(fs_init(1)+fs_init(2))*(1+gi)+Num/Fs2*fs_init(3)*(0.5*gi));  ones(1,Nc)*(2*Num/Fs2*(fs_init(1)+fs_init(2)+fs_init(3))*(1+gi)+Num/Fs2*fs_init(4)*(0.5*gi))]; %% three bands  
index=round(sub_band*([1:1:Nc]-1)+1+offset_4); 
Pm=H_freq(index) ;                                                        %%  sub channel gain in DMT
p_sq=(sum(abs(Pm').^2))'/Nc;

SNR_MFB=p_sq*Epsilon_x./(Sigma_r+PARr*p_sq*Epsilon_x*2^(-2*b0).*fs_init.^(2*bs));
Gama_DMT=p_sq./((prod(abs(Pm').^2))'.^(1/Nc)) %%  For matrices, prod(X) is a matrix the same size as X containing the cumulative products over each column. 
% ((prod(abs(Pm').^2))'.^(1/Nc))
C_DMT=log2(SNR_MFB./(Gama_DMT*Gap))     
Used_band=sum(fs_init*(1+ai)*(1+gi));
Bitrate=sum(fs_init.*C_DMT)/1e3  %Gbps


%% Optimization 
%Bitrate=sum(fs_init.*(log2((p_sq*Epsilon_x./(Sigma_r+PARr*p_sq*Epsilon_x*2^(-2*b0).*fs_init.^(2*bs)))./((p_sq./((prod(abs((H_freq(round(sub_band*([1:1:Nc]-1)+1+[ones(1,Nc)*Num/Fs2*fs_init(1)*gi; ones(1,Nc)*(2*Num/Fs2*fs_init(1)*(1+gi)+Num/Fs2*fs_init(2)*gi)])))').^2))'.^(1/Nc)))*Gap))))/1e3

 Fun=@(fs_init)1e6/(sum(fs_init.*(log2((p_sq*Epsilon_x./(Sigma_r+PARr*p_sq*Epsilon_x*2^(-2*b0).*fs_init.^(2*bs)))./((p_sq./((prod(abs((H_freq(round(sub_band*([1:1:Nc]-1)+1+offset_4)))').^2))'.^(1/Nc)))*Gap)))))
% 
 A=[(1+ai)*(1+gi) (1+ai)*(1+gi) (1+ai)*(1+gi) (1+ai)*(1+gi)];
 Aeq=[];
 Beq=[];
 LB=[1;1]; %% FS_min> 1M
 UB=[];
  %options = optimset('Algorithm','active-set'); 
options = optimset('Algorithm','sqp');
[X Bitrate]= FMINCON(Fun,[20;50;100;30],A,Band_total,Aeq,Beq,[1,1,1,1],UB,[],options);
X
1e3./Bitrate

% =======================The following generate the bits distribution vs frequency=====================
offset_opm=[ones(1,Nc)*Num/Fs2*X(1)*0.5*gi;  ones(1,Nc)*(2*Num/Fs2*X(1)*(1+gi)+Num/Fs2*X(2)*0.5*gi);  ones(1,Nc)*(2*Num/Fs2*(X(1)+X(2))*(1+gi)+Num/Fs2*X(3)*(0.5*gi)); ones(1,Nc)*(2*Num/Fs2*(X(1)+X(2)+X(3))*(1+gi)+Num/Fs2*X(4)*(0.5*gi))];
sub_band=2*Num*X/Fs2/Nc; 
index_opm=round(sub_band*([1:1:Nc]-1)+1+offset_opm); 
Pm=abs(H_freq(index_opm)).^2;
sub_band=2*X/Fs2/Nc;    
%when N changes, should modefy dimension matrix of X_matrix
X_matrix(1,1:32)=X(1);
X_matrix(2,1:32)=X(2);
X_matrix(3,1:32)=X(3);
X_matrix(4,1:32)=X(4);
SNR_wADC=Epsilon_x*Pm./(Sigma_r+PARr*Pm*Epsilon_x*2^(-2*b0).*X_matrix.^(2*bs));
%SNR_wADC=Epsilon_x./(Sigma_r./Pm+PARr*Epsilon_x*2^(-2*b0)*X_matrix.^(2*bs));
SNR_woADC=Epsilon_x*Pm./(Sigma_r);
Bits_wADC=log2(SNR_wADC/Gap);     %bits/Hz/s
Bits_woADC=log2(SNR_woADC/Gap);     %bits/Hz/s
figure(2)
plot([400/4096*index_opm(1,1)+0.01,400/4096*index_opm(1,:),400/4096*index_opm(1,32)+0.01,...
    400/4096*index_opm(2,1)-0.001,400/4096*index_opm(2,:),400/4096*index_opm(2,32)+0.01,...
    400/4096*index_opm(3,1)-0.001,400/4096*index_opm(3,:),400/4096*index_opm(3,32)+0.01,...
    400/4096*index_opm(4,1)-0.001,400/4096*index_opm(4,:),400/4096*index_opm(4,32)+0.01],...
    [0,Bits_woADC(1,:),0,...
    0,Bits_woADC(2,:),0,...
    0,Bits_woADC(3,:),0,...
     0,Bits_woADC(4,:),0],'--b','LineWidth',3);
hold on;
plot([400/4096*index_opm(1,1)+0.01,400/4096*index_opm(1,:),400/4096*index_opm(1,32)+0.01,...
    400/4096*index_opm(2,1)-0.001,400/4096*index_opm(2,:),400/4096*index_opm(2,32)+0.01,...
    400/4096*index_opm(3,1)-0.001,400/4096*index_opm(3,:),400/4096*index_opm(3,32)+0.01,...
    400/4096*index_opm(4,1)-0.001,400/4096*index_opm(4,:),400/4096*index_opm(4,32)+0.01],...
    [0,Bits_wADC(1,:),0,0,Bits_wADC(2,:),0,0,Bits_wADC(3,:),0,...
    0,Bits_wADC(4,:),0],'k','LineWidth',3);
title('4 Bands','fontsize',14)
xlabel('Frequency(MHz)','fontsize',14);
ylabel('Bits','fontsize',14);
legend('Without ADC noise','With ADC noise');
set(gca,'fontsize',14)
hold on


%============================contour optimization======================
ff_1=(1:3:380)/1.2;
ff_2=(1:3:380)/1.2;
for i=1:length(ff_1)
for j=1:length(ff_2)
if((ff_2(i)*1.2+ff_1(j)*1.2)>400)
    Bitrate_contour(i,j)=0;
else
offset=[ones(1,Nc)*Num/Fs2*ff_1(j)*gi; ones(1,Nc)*(2*Num/Fs2*ff_1(j)*(1+gi)+Num/Fs2*ff_2(i)*gi)];
sub_band=2*Num*[ff_1(j);ff_2(i)]/Fs2/Nc;
index=round(sub_band*([1:1:Nc]-1)+1+offset); 
Pm=H_freq(index) ;                                                        %%  sub channel gain in DMT
p_sq=(sum(abs(Pm').^2))'/Nc;
Bitrate_contour(j,i)=sum([ff_1(j);ff_2(i)].*(log2((p_sq*Epsilon_x./(Sigma_r+PARr*p_sq*Epsilon_x*2^(-2*b0).*[ff_1(j);ff_2(i)].^(2*bs)))./((p_sq./((prod(abs((H_freq(round(sub_band*([1:1:Nc]-1)+1+offset)))').^2))'.^(1/Nc)))*Gap))));
end
end
end
% figure(3)
% [c,h]=contour(ff_2,ff_1,Bitrate_contour,15);
% clabel(c,h,'manual');
% xlabel('Band2(MHz)','fontsize',14);
% ylabel('Band1(MHz)','fontsize',14);

