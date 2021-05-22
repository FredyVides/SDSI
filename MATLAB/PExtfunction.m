function y=PExtfunction(x,N)
	H=@(x,a)(x>a);
	y=0;
	for k=0:(N-1)
		y=y+min(x-k,1-x+k).*(H(x,k)-H(x,k+1));
	end
end
