%GOPH 547-Gravity and Magnetics 
%Safian Omar Qureshi
%ID: 10086638
%TA: Ye Sun
 
%% Question 1a
function [U]=grav_pot_point(x,xm,m,G) %defining gravity potential function

r_norm=norm(x-xm);                    %simply using formula 
U=(G*m)./r_norm;

end


%these will be used in other parts of the code 

