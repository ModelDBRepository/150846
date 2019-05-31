% This function computes the sigmoid transfer function used by Model_GridCell.m and by PlastRuleDemo.m
% Luisa Castro, FCUP
% luisa.castro@fc.up.pt

function[F]=SGtransf(I,a,b)
    
    % Sigmoid transfer function with parameters:
    %   a - gives the width of the slope (high values make the slope smoother)
    %   b - gives the center (inflexion) point of the slope (input units)
    
    F=1./(1+exp(-(I-b)/a));


