% exactfun.m
% Peter Ferrero, Oregon State University, 4/15/2018, MTH 553, HW 2
% A function to calculate the exact solution of the 2D Poisson's equation

function uExact = exactfun(x,y)

uExact = sin(pi*x)*sin(pi*y);

end