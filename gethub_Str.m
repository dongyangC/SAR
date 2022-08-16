%% Strimap_SAR Mode
clear all;clc;close all;
%% 参数
C = 2.99792458e+8;        
fc = 5.400012e9;             
lambda = C/fc;             
Va = 7568.945978;  
PRF =1410.025390;               
Tr =20e-6;                
Br = 60e6;                
Fr = 72e6;                
Kr = Br/Tr;              
R0 = 815900.637233;              
D  = 15;                   
theta_rc=0/180.*pi;       
Ls = R0*lambda/D;                                                                                                             
Tsar = Ls/Va;                                                                                                 
%% 
Dx = D/2;                
Dr = C/Br/2;                                                                  
%% 
Na=ceil((2*Ls/Va)*PRF); 
Na=2048;
ta=(-Na/2:Na/2-1)/PRF; 
%% 
r0=1000;                            
dt=1/Fr;                            
Rmin=sqrt((R0-r0).^2-(Ls./2).^2);         
Rmax=sqrt((R0+r0).^2+(Ls/2).^2); 
Nr=ceil((2*(Rmax-Rmin)/C/dt+Tr/dt)./2);
Nr=2^nextpow2(2*Nr);                
tr = 2*R0/C+(-Nr/2:Nr/2-1)/Fr;
%%
Ptarget=[  
          0,R0,1;
        ];                 
Ntarget =3;
T=Ptarget; 
Raw = zeros(Na,Nr);
%%
sita_r_c = (2*pi)/180;
for k=1:1
    sigma=T(k,3);
    Dslow=ta*Va-T(k,1);   
    Rr=sqrt(Dslow.^2+T(k,2).^2);
    tau=2*Rr/C;
    Dfast=ones(Na,1)*tr-tau'*ones(1,Nr);
    phase=pi*Kr*Dfast.^2-(4*pi/lambda)*(Rr'*ones(1,Nr));
    Raw=Raw+sigma*exp(j*phase).*(abs(Dfast)<=Tr./2).*((abs(Dslow)<=Ls./2)'*ones(1,Nr));
end
%% 
A_ph_3 = [exp(j*120*pi/180),exp(j*150*pi/180),exp(j*175*pi/180)];
A_NN =3; 
aaaa=[1,1,1];
Raw1=zeros(Na,Nr);
Raw2=zeros(Na,Nr);
 for i=1:A_NN
    Raw(:,i:A_NN:end)=Raw(:,i:A_NN:end).*A_ph_3(1,i);
 end
A_ph_3 = [exp(j*180*pi/180),exp(j*150*pi/180),exp(j*120*pi/180)];
R_NN =3;
for i=1:R_NN
  Raw(i:R_NN:end,:)= Raw(i:R_NN:end,:).*A_ph_3(1,i);
end
%% 
fr=Fr/Nr*(-Nr/2:Nr/2-1);
fa=PRF/Na*(-Na/2:Na/2-1); 
alpha = 1;
beta = (1 - (fa*lambda/2/Va).^2).^0.5;   
a = 1./beta - 1;     
R = R0./beta;         
a_scl = a + (1-alpha).*(1+a)./alpha;    %距离变标系数
k_inv = -1./Kr - (2.*lambda.*R0.*(beta.^2-1))./(C^2.*beta.^3);
Km = 1./k_inv;       
data = ftx(Raw);
H1 = exp(-j*pi*((Km.*a_scl)'*ones(1,Nr).*((ones(Na,1)*tr- (2.*R./C)'*ones(1,Nr)).^2)));
data = data .* H1;
data = fty(data);
H2 = exp(-j*pi*((1./(Km.*(1+a_scl)))'*ones(1,Nr)).*(ones(Na,1)*fr.^2)) .* exp(j*4*pi*R0/C.*((ones(Na,1)*fr) .* (a'*ones(1,Nr))));
data = data .* H2;
data = ifty(data);
R_0 = tr/2*C;
dphi = 4*pi*((Km.*a_scl.*(1+a).^2 ./ (C^2.*(1+a_scl)))'*ones(1,Nr)).*(ones(Na,1)*((R_0-R0).^2));
H3 = exp(j*dphi);
data = data .* H3;
H4 = exp(j*4*pi/lambda*(ones(Na,1)*R_0).*((beta-1)'*ones(1,Nr)));
data = data .* H4;
clear H4;
data = iftx(data);
%%
tt=20*log10(abs(data)./max(max(abs(data))));
figure
imagesc(tt,[-75,0]);
colormap(jet);