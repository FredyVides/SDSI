tic,Udata=importdata('../DataSets/GFUdata.csv');toc

[d,r,~,U0,U1]=deg(3,1,Udata(:,1:400),9e-7);
size(U0)


disp('=========================================================');
disp('Computing sparse identification with SINDy:');
tic,A=SINDy(U0,U1,500,6e-3);toc
disp('=========================================================');
disp('=========================================================');
disp('Computing sparse identification with SpSolver:');
tic,A0=SpSolver(U0,U1,r,7e-3);toc
disp('=========================================================');

subplot(221),spy(A)
subplot(222),spy(A0)

tic,                          
u1=U0(:,1);
u2=u1;
v0=zeros(r-1,1);
v0(1)=1;
M=size(Udata,2); 
Ak=A0;
for k=1:(M-1), u1=[u1 U0*(Ak*v0)];Ak=A0*Ak;end
csvwrite('../DataSets/PGFUdata0.csv',u1);
  
Ak=A;
for k=1:(M-1), u2=[u2 U0*(Ak*v0)];Ak=A*Ak;end
csvwrite('../DataSets/PGFUdata1.csv',u2);
toc

tic,l=eig(A);toc
tic,l0=eig(full(A0));toc

subplot(223),plot(real(l),imag(l),'r.'),axis square;axis tight;
subplot(224),plot(real(l0),imag(l0),'r.'),axis square;axis tight;
