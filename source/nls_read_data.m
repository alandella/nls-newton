clc
close all
clearvars

%% scalar equation

x0 = 2.1;

options = optimoptions('fsolve','display','iter');

xf = fsolve(@test,x0,options)

Fx = test(xf)

%% vector system

dim = 3;

x0_vec = ones(dim,1);
 
%options = optimoptions('fsolve','display','iter');
 
xf = fsolve(@(x) test_vec(x,dim),x0_vec,options)

Fx = norm(test_vec(xf,dim))

%% function declaration

function f = test(x)

    f = 5*x + log(abs(x)) - 10000;

end

function f = test_vec(x,dim)

    f = zeros(dim,1);
    
	f(1) = 3*x(1) - cos(x(2)*x(3)) - 3/2;
	f(2) = 4*x(1).^2 - 625*x(2).^2 + 2*x(3) - 1;
	f(3) = 20*x(3) + exp(-x(1).*x(2)) + 9;
    
end

