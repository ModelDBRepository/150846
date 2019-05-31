function [ r ] = PC_2DGaus(x,mu,noise)
%PlaceField_2DGaussian Calculate the firing rate of a place cell using a generalized 2D Gaussian model
%   Place field parameters
%   px   : x position of the animal 
%   py   : y position of the animal
%   x0  : x center of place field
%   y0  : y center of place field
%   sx  : scale in x (std_x)
%   sy  : scale in y (std_y)
%   rho : angle of ellipse cross-section with xx axis
%   noise : value of noise between 0 and 1 
%(note that x =[px py] and mu =[xo yo sx sy rho ]
% Authors: Paulo Aguiar and Luisa Castro

px=x(:,1);
py=x(:,2);
x0=mu(:,1);
y0=mu(:,2);
sx=mu(:,3);
sy=mu(:,4);
rho=mu(:,5);

sx2 = 2.0 * sx .* sx;
sy2 = 2.0 * sy .* sy;

a =  cos(rho).^2   ./ (sx2)     + sin(rho).^2   ./ sy2;
b = -sin(2.0*rho) ./ (2.0*sx2) + sin(2.0*rho) ./ (2.0*sy2);
c =  sin(rho).^2   ./ (sx2)     + cos(rho).^2   ./ sy2;

r = exp( -( a.*(px-x0).^2 + 2.0*b.*(px-x0).*(py-y0) + c.*(py-y0).^2 ) );

r = r.*(1+noise*randn(length(r),1));

r = max(r,0);
end

