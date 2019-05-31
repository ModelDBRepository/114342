
%  FD weights  case 2 (using all grid points)



function [DIFF,x]=dmc(N,M,b)
sf=1/cos(pi/N);

% xi-----point at which the approximation is to be correct-----------
% 

for II=0:N,

	xi=-sf*b*cos(II*pi/N)+b;

% M----highest order of derivative
%M=3;

% x-----Chebyshev exteme nodes on the interval [-d,2*b+d]
x=zeros(N+1,1);
for i=0:N,
	x(i+1)=-sf*b*cos(pi*i/N)+b;
end

c=zeros(N+1,M+1);

c(0+1,0+1)=1.0;
c1=1.0;
c4=x(0+1)-xi;

for I=1:N,
	MN=min([I,M]);
	c2=1.0;
	c5=c4;
	c4=x(I+1)-xi;
	
	for J=0:I-1,
		c3=x(I+1)-x(J+1);
		c2=c2*c3;
		
		for K=MN:-1:1,
			c(I+1,K+1)=c1*(K*c(I,K)-c5*c(I,K+1))/c2;
		end

		c(I+1,1)=-c1*c5*c(I,1)/c2;
		
		for K=MN:-1:1,
			c(J+1,K+1)=(c4*c(J+1,K+1)-K*c(J+1,K))/c3;
		end
		
		c(J+1,1)=c4*c(J+1,1)/c3;
	end
  c1=c2;
end

DIFF(II+1,:)=c(:,M+1)';

end

%xi
%x

