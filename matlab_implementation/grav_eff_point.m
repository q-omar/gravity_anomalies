%GOPH 547-Gravity and Magnetics (W2017)
%Safian Omar Qureshi
%ID: 10086638
%TA: Ye Sun

%% Question 1b
function [gz]=grav_eff_point(x,xm,m,G) %defining gravity effect function

r_norm=norm(x-xm);                     %simply using formula
r3=r_norm^3;                               
dz=(x(3)-xm(3));
gz=(G*m*dz)/r3;

end

%these will be used in other parts of the code 
