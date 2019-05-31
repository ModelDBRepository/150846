% This code generates a trajectory for the rat
% within the square maze, for the main model Model_GridCell.m
% Luisa Castro, FCUP
% luisa.castro@fc.up.pt

mean_v=0.08;						    % m/s, rat's mean velocity
v=mean_v+0.0*randn(1,bins);  			% m/s, vector of rat's velocity for the rat walk
h=1000./(v*dt);       				    % /m, adjusted in order to have a mean velocity of v (1/h is the distance travelled in each time step)
x=zeros(bins,2);     					% initialization of positions vector
x(1,:)=[x1 y1];      					% m, initial point of rat walk, defined in main model Model_GridCell.m
alf=zeros(1,bins);   					% radians, stores the effective directions of movement

direc=pi/120; 						    % tuned according to dt and v

for t=2:T/dt
    alf(t)=alf(t-1)+direc*randn;
    x(t,:)=x(t-1,:)+[cos(alf(t))/h(t) sin(alf(t))/h(t)];
    while (x(t,1)  >  side  ||   x(t,2)   >  side ||  x(t,1)   <   0   ||   x(t,2)  <   0 )
        alf(t)=alf(t-1)+(pi/4)*randn;		
        x(t,:)=x(t-1,:)+[cos(alf(t))/h(t) sin(alf(t))/h(t)]; % the virtual rat makes an aprox 90deg turn if he hits the wall
    end  
end

figure			%Fig. 2b
plot(x(:,1),x(:,2))
title('Trajectory of the virtual rat')
xlim([0 side]); xlabel('m'); set(gca,'XTick',[0 1])
ylim([0 side]); ylabel('m'); set(gca,'YTick',[0 1]); axis square

% Computing path's length = travelled distance and mean velocity 
dx=diff(x(:,1));
dy=diff(x(:,2));
path_m=sum(sqrt(dx.*dx+dy.*dy));   		% m
velocity_ms=path_m/(T/1000) 			% m/s
