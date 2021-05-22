function Lags=LagEstimate(signal,threshold)
	[r,lags] = xcorr(signal);
	r=abs(r);
	r=r/max(r);	
	d0=diff(r(1:(end-1)));
	d1=diff(r(2:end));
	f=find(d0.*d1<0);
	g=find(abs(r(f+1))>threshold);
	Lags=lags(f(g)+1);
	f=find(Lags>0);
	Lags=Lags(f);
end
