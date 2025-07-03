function y = Ackley(x) 
% Ackley function
global NFE
NFE = NFE + 1;


d = numel(x);
a = 20; b = 0.2; c = 2*pi;
s1 = 0; s2 = 0;
for i=1:d
   s1 = s1+x(i)^2;
   s2 = s2+cos(c*x(i));
end
y = -a*exp(-b*sqrt(1/d*s1))-exp(1/d*s2)+a+exp(1);
end