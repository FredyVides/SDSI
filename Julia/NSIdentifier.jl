using DataFrames,CSV,LinearAlgebra,BenchmarkTools,Plots,DelimitedFiles

include("spsolver.jl")

@elapsed df=CSV.read("../DataSets/GFUdata.csv",DataFrame;header=0) |> Matrix{Float64}

rk=380;

df0=df[:,1:rk];
df1=df[:,2:(rk+1)];

@elapsed A=df0\df1

@elapsed A0=spsolver(df0,df1,rk,5e-3,5e-3)

u0=df0[:,1];                            
u1=u0;
u0=df0\u0;  
M=size(df,2); 
Ak=A0;
for k=1:(M-1)
    global u1=[u1 df0*(Ak*u0)];
    global Ak=A0*Ak;
end
writedlm("../DataSets/PGFUdata0.csv",u1,',');
u2=df0*u0;  
Ak=A;
for k=1:(M-1) 
    global u2=[u2 df0*(Ak*u0)];
    global Ak=A*Ak;
end
writedlm("../DataSets/PGFUdata1.csv",u2,',');

@elapsed l=eigvals(A)

@elapsed l0=eigvals(A0)

scatter(real(l),imag(l))

scatter(real(l0),imag(l0))
