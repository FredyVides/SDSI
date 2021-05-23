% [H,H0,H00,w0,c,c0,c00]=NLSESpModelID(dt,data,difforder,L,tol,Ls,delta)
%%
% Approximate sparse model identification
% 
% Example:
%
% Import and preprocess data (data generated with NLSchrodigner1DRK.m):
% nwdata=csvread('../DataSets/NoisyNLSEqData.csv');
% rwdata=csvread('../DataSets/NLSEqData.csv');
% tol=8e-4;Threshold= 2e-2;noise=1e-6;
% [M,N]=size(nwdata);
% x=nwdata(:,1)';
% Wdata=rwdata(:,2:N);
% wdata=nwdata(2:(M-1),2:N);
% dt=.01;
% Compute sparse model identification
% L=200;S=35;
% tic,[H,H0,H00,w0,c,c0,c00]=NLSESpModelID(2*i*dt,wdata,4,S,tol,L,Threshold);toc
% subplot(211),spy(H),subplot(212),stem(abs(c))
% Compute and visualize predictions
% [w,wdata1]=NLSEPredictor(dt,2,w0,x,H,c(4:203),1000);
% Estimate and visualize prediction errors:
% Error0=max(abs(Wdata(:,S:(S+500))-w));
% Error1=max(abs(abs(Wdata(:,S:(S+500)))-abs(w))); 
% figure;
% subplot(211),semilogy(Error0),axis tight,grid on;
% subplot(212),semilogy(Error1),axis tight,grid on;
%
% Author: Fredy Vides <fredy@HPCLAB>
% Created: 2021-01-06
%%
%%
function [H,H0,H00,w0,c,c0,c00]=NLSESpModelID(dt,data,difforder,L,tol,Ls,delta)
M=size(data,1);
data=data(:,1:L);
dtdata=DtData(data,dt,difforder);
dtdata=dtdata(:);
E=@(j,n)double((1:n)==j);
Z=zeros(1,M);
Sdata=data(:);
for k=2:(difforder-2)
    S=toeplitz(E(k,M),Z);
    sdata=S*data;
    Sdata=[Sdata sdata(:)];
    sdata=S'*data;
    Sdata=[Sdata sdata(:)];
end
for k=1:Ls
	sdata=(abs(data).^k).*data;
	Sdata=[Sdata sdata(:)];
end
Sdata=Sdata;
dtdata=dtdata;
disp('=========================================================');
disp('Computing sparse identification with SpSolver:');
tic,c=SpSolver(Sdata,dtdata,Ls,tol,delta);toc
disp('=========================================================');
tau = 1e-5; mu=1/2; MaxIt = 1e5; tol= 1e-12;sigma=2e-2;
disp('=========================================================');
disp('Computing sparse identification with SINDy:');
tic,c0=SINDy(Sdata,dtdata,500,sigma);toc
disp('=========================================================');
disp('=========================================================');
disp('Computing sparse identification with Douglas Rachford:');
tic, c00=DouglasRachford(Sdata,dtdata,sigma,tau,mu,MaxIt,tol);toc
disp('=========================================================');
H=c(1)*toeplitz(E(1,M));
H0=c0(1)*toeplitz(E(1,M));
H00=c00(1)*toeplitz(E(1,M));
Lk=max(2,(2*(difforder-3)-1));
for k=2:2:Lk
    S=toeplitz(E(k,M),Z);
    H=H+c(k)*S+c(k+1)*S';
    H0=H0+c0(k)*S+c0(k+1)*S';
    H00=H00+c00(k)*S+c00(k+1)*S';
end
w0=data(:,L);
end
