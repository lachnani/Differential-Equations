function [x] = rk4system(eqnum,h,t,x,n)
%RK4system Solves system of ODEs unsing RK method
%   Takes in number if equations, step size, initial t, function array, 
%   and the number of steps.
y = zeros(eqnum,1);
k = zeros(eqnum,4);
for j = 1:n
    k(:,1) = xpsystem(t,x,k(:,1));
    for i = 1:eqnum
        y(i) = x(i) + h*k(i,1)/2;
    end
    k(:,2) = xpsystem(t+(h/2),y,k(:,2));
    for i = 1:eqnum
        y(i) = x(i) + h*k(i,2)/2;
    end
    k(:,3) = xpsystem(t+(h/2),y,k(:,3));
    for i = 1:eqnum
        y(i) = x(i) + h*k(i,3);
    end
    k(:,4) = xpsystem(t+h,y,k(:,4));
    for i = 1:eqnum
        x(i) = x(i) + h*(k(i,1)+2*k(i,2)+2*k(i,3)+k(i,4))/6;
    end
    t = t+h;
end
end

