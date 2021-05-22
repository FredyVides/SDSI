function y=LTIPredictor(A,y0,N)
	y0=y0(:);
	y=y0.';
	y0=flipud(y0);
	L=length(y);
	for k=1:N
		z0=A*y0;
		y0=circshift(y0,1);
		y0(1)=z0;
		y=[y z0];
	end
end
