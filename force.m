function [F] = force(A,B)
global sigma k0
r=((A(1)-B(1))^2+(A(2)-B(2))^2+(A(3)-B(3))^2)^0.5;
e=[(A(1)-B(1))/r,(A(2)-B(2))/r,(A(3)-B(3))/r];
f=k0*(r-0.5*sigma);
F=f*e;
end