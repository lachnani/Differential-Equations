function [x] = eulermess(a,b,xa,n)
%EULERMESS computes taylor approximation of order v for an ODE
%   Must be paired with feuler=f(t,x) 
h=(b-a)/n;
x=xa;
%organize the output
% fprintf('\n')
% disp('                     Euler Method')
% disp('___________________________________________________________')
% disp('   ti      f(ti,xi)        xi         Exact         Error')
% disp('___________________________________________________________')
% fprintf('\n')
% fprintf('%6.2f       ---    %12.6f %12.6f      %4.2f\n',a,x,x,0)
%main loop
for i=1:n
    t=a+(i-1)*h;
    m=feuler(t,x);
    x = x+h*m(1)+ (h^2)*m(2)/2 + (h^3)*m(3)/6 +(h^4)*m(4)/20;
    %Display results
    t=t+h;
    y=feulersol(t);
    err=abs(y-x);
%     fprintf('%6.2f %12.6f %12.6f %12.6f      %8.2e\n',t,m(1),x,y,err)
end    
end