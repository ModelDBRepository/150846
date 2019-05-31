%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Code to produce main figures/results from:                                    %
%                                                                               %
% " A feedforward model for the formation of a grid field where spatial         %
%   information is provided solely from place cells "                           %
% Biological Cybernetics, 2014                                                  %
%                                                                               %
% Authors:                                                                      %
%   Luisa Castro and Paulo Aguiar                                               %
%   Faculdade de Ciencias da Universidade do Porto, Portugal                    %
%   contact email: pauloaguiar@fc.up.pt                                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Simulation running, please wait. Press ctrl+c to end... \n')
close all
clear all

%% Boolean and noise parameters to choose

plot_online=0;      % use 1 for online plot during the simulation (slow) ; 0 for no online plot
rand_PC_centers=0;  % use 0 for homogeneous scatter over a square lattice and 1 for random locations of PC's centers
anisot_PC=0;        % use 0 for isotropic place fields or 1 for different shape place fields
use_rate=0;         % use 1/0 to evaluate grid score from rate grid map/total input
sim_rat_walk=1;	    % use 1 for new rat walk ; 0 for uploading rat walk (as file x.mat)
x_haft=0;           % use 1 if using Hafting's 2008 rat path

rate_noise=0.0;     % use 0 for no noise in PC firing rates or above (below 1) for additive Gaussian noise

%% Time of the simulation
Tmin=20;      	    % min, time duration of the simulation (adapt to 29.882 if using imported rat trajectory)
T=Tmin*60*1000;     % ms, time duration of the simulation

dt=1;               % ms, step size (change to 20ms if using Hafting path)
side=1;             % m, square maze side length

times_sim=0: dt:T-dt;               % ms, simulation steps [ms]
bins=length(times_sim);             % number of steps in the simulation

%% Rat Walk Trajectory

if ~sim_rat_walk
    if x_haft
        %The file x_haft3 has the following rescaling of Hafting et al 2008 rat positions (pos_x,pos_y):
        %x_haft=(90+[pos_x pos_y])./100;
        load x_haft3
        x=x_haft3;
        figure			   % Fig. 2a
        plot(x(:,1),x(:,2)); title('Example of a real rat trajectory [Fig. 2a]')
        side = 1.8;
        xlim([0 side]); xlabel('m'); set(gca,'XTick',[0 1])
        ylim([0 side]); ylabel('m'); set(gca,'YTick',[0 1]); axis square
    else
        load x          % upload the x rat trajectory
    end
    
else
    x1=0.5;             % m, 1st coordinate of the initial position of the rat
    y1=0.5;             % m, 2nd coordinate of the initial position of the rat
    rat_walk            % generate a new rat trajectory
end

%% Place Cells Parameters

Sigma = 0.05;                                                      % m, place field (pf) width to produce dorsal fields with approx. 0.12m radius
dx_mu=1/100;                                                       % place cells density
add_bord=0.4;                                                      % m, additional border for distribute place cells centers
big_side=side+2*add_bord;                                          % m, side of big maze = available maze for the rat + additional border
aux1mu=0-add_bord+dx_mu/2:dx_mu:side+add_bord-dx_mu/2;             % m, 1st coord of pfs centers - for convenience, no cells on borders
aux2mu=0-add_bord+dx_mu/2:dx_mu:side+add_bord-dx_mu/2;             % m, 2nd coord of pfs centers - for convenience, no cells on borders
mux=repmat(aux1mu,1,length(aux2mu));
muy=sort(repmat(aux2mu,1,length(aux1mu)));
mu_cent=[mux' muy'];                                               % m, place fields 2d centers in vector Nx2

N_in=length(find(and(and(mu_cent(:,1)>=0,mu_cent(:,1)<=side),and(mu_cent(:,2)>=0,mu_cent(:,2)<=side))));      % count cells with pf center inside the maze
N=length(mu_cent);                                                                             % count all place cells

if rand_PC_centers
    figure;plot(mux,muy,'.');axis square
    mux_rand=dx_mu.*rand(1,N)+mux-dx_mu/2;
    muy_rand=dx_mu.*rand(1,N)+muy-dx_mu/2;
    mu_rand=[mux_rand' muy_rand'];
    mu_cent=mu_rand;
end

figure  %Fig. 7a or 7b
plot(mu_cent(:,1),mu_cent(:,2),'.');xlim([0 0.2]);ylim([0 0.2]); axis square   
title( 'place fields spatial distribution' );

if ~anisot_PC
    ro=zeros(N,1);
    Sigmax=Sigma*ones(N,1);
    Sigmay=Sigma*ones(N,1);
else
    prod=Sigma^2*(1+0.1*randn(N,1));   % normal distribution for area of fields with mean Sigma^2
    rac=1+0.1*randn(N,1);              % normal distribution for axis ratio with mean 1
    Sigmax=sqrt(prod.*rac);
    Sigmay=prod./sqrt(prod.*rac);
    ro=pi*rand(N,1);                   % angle in radians for the ellipse axis rotation
end

mu=[mu_cent Sigmax Sigmay ro];

scale=side*side/(N_in*2*pi*Sigma*Sigma);     % scaling parameter for pc's weights: for unitary weights, input to GC becomes 1

%% Grid Cell Parameters

tau_r=30;    % [ms] time constant for the output firing rates

%% Weight and their Leaning Rule Parameters

r_th=0.98;     % [normalized rate units] threshold value for updating weights rule

gain=2.5;
sigmatil=gain*Sigma;			  % width for gain function
mug=[mu(:,1:2)    gain*mu(:,3)     gain*mu(:,4)        mu(:,5)];

w_b=0.97*r_th*scale;              % weight value for baseline contributing cells
w_s=0.99*r_th*scale;              % weight value for R3
wmax=1.5*r_th*scale;              % maximum weight for each synapse, for R1 and R3 intersections
wmin=0;					          % minimum weight for each synapse, for R2

MW = w_b*ones(N,1);               % initialization of weights matrix

% Plasticity parameters for dorsal zone :=  Grid Cells fields: radius=0.10 m and spacing G=0.35 m
phi1=0.07;                        % m, R1 radius
phi2=0.26;      		          % m, inner radius for R3
phi3=0.34;              		  % m, outer radius for R3

%%SUBSTITUIR OS SEGUINTES VALORES USANDO A FUNCAO INT_FUNC_PC?
q1=exp(-phi1^2/(2*(sigmatil)^2));
q2=exp(-phi2^2/(2*(sigmatil)^2));
q3=exp(-phi3^2/(2*(sigmatil)^2));

%% Cartoon plot for plasticity weight level for each value of the gain function (Fig.3)
resol=0.001;
x_plot=0:resol:1;
frate= PC_2DGaus([x_plot'  0*x_plot'],[0.5 0 Sigma Sigma 0],rate_noise);              % firing rate profile
gfrate=PC_2DGaus([x_plot'  0*x_plot'],[0.5 0 gain*Sigma gain*Sigma 0],rate_noise);   % gain function of firing rate profile
w_plot=w_b*ones(1,1+1/resol);
w_plot(gfrate>=q1)=wmax;
w_plot(and(gfrate<q1,gfrate>=q2))=wmin;
w_plot(and(gfrate<q2,gfrate>=q3))=w_s;
figure
[AX,H1,H2]=plotyy(x_plot,gfrate,x_plot,w_plot);
set(get(AX(1),'Ylabel'),'String','Gain Function');     set(AX(1),'YTick',[0 1])
set(get(AX(2),'Ylabel'),'String','Synapse Weight');
set(AX(2),'YTick',[wmin w_b w_s wmax]);
set(AX(2),'YTicklabel',{'wmin','w_b','', 'wmax'});
set(H1,'LineWidth',3);                    			  set(H2,'Marker','.')
set(AX,'XTick', [0 1]);					       xlabel('Distance [m]')

%% VARIABLES STORAGE AND INITIALIZATION

firrat_pc=zeros(N,1);               % firing rates of place cells in each time
firrat_gr=zeros(1,bins);            % firing rate of output cell in each time

%% MAIN CYCLE
tic

if plot_online
    h1=figure;
end

nfields=0;                            % keep track of the number of grid fields created for the GC
no_plast_zones=0;                     % turns 1 if no more plasticity zones exist in the available maze
acetil =0.2*r_th*scale*ones(N,1);     % weight increase due to acetylcoline effect

for i=2:T/dt
    
    if nfields>=2
        acetil=0*acetil;              % after the second node is formed, abolishes the acetilcoline effect
    end
    
    t=(i-1)*dt;
    
    firrat_pc(:)=PC_2DGaus(x(i,:),mu,rate_noise);
    
    I= (min(MW+acetil.*(MW>0),wmax))'*firrat_pc(:);
    
    firrat_gr(i)=firrat_gr(i-1)+(dt/tau_r)*(-firrat_gr(i-1)+SGtransf(I,0.12,0.5));     % Computes firing rate value for Grid cell
    
    % Weight Modification
    
    if firrat_gr(i)>=r_th
        gfr=PC_2DGaus(x(i,:),mug,rate_noise);             % Gain of PCs firing rate
        
        MW(gfr>=q1)=wmax; 				                	                 % R1 - central excitatory disc
        MW(    and(     gfr<q1  ,   gfr>=q2))=wmin;                   		 % R2 - inhibitory ring
        MW(    and(     gfr>=q3 ,   and(MW(:)<wmax,MW(:)>=w_s)  ))=wmax;	 % intersection diamond zones
        MW(    and(     gfr>=q3 ,   and(MW(:)<=w_s,MW(:)>0)  ))=w_s;         % R3 - excitatory ring
        
        nfields=nfields+1
        
        if plot_online
            figure (h1)
            
            subplot(1,2,1)
            MW1f_xy=reshape(MW(:),sqrt(N),sqrt(N));
            imagesc(MW1f_xy)
            set(gca,'YDir','normal')
            set(gca,'XTick',[]); set(gca,'YTick',[]); xlabel('m');   ylabel('m')
            axis square
            title({'Place Cell Weights';'(available maze)'})
            xlim([(sqrt(N)-sqrt(N_in))/2 sqrt(N_in)+(sqrt(N)-sqrt(N_in))/2])
            ylim([(sqrt(N)-sqrt(N_in))/2 sqrt(N_in)+(sqrt(N)-sqrt(N_in))/2])
            subplot(1,2,2)
            x1 =0+dx_mu/2:dx_mu:side-dx_mu/2; x2 =0+dx_mu/2:dx_mu:side-dx_mu/2;
            ExampInp=zeros(length(x1),length(x2));
            for idpc=1:N
                [X1,X2] = meshgrid(x1,x2);
                Fplotw =MW(idpc)*PC_2DGaus([X2(:) X1(:)],mu(idpc,:),rate_noise);
                Fplotw = reshape(Fplotw,length(x2),length(x1));
                ExampInp=ExampInp+Fplotw;
            end
            imagesc(ExampInp)
            set(gca,'XTick',[]); set(gca,'YTick',[]); xlabel('m');   ylabel('m')
            set(gca,'YDir','normal') ;
            axis square
            title({'Input(\Sigma w*r) to GC'; '(available maze)'})
        end
        
    end
    
    
    %Stop if no more Input>r_th exists (evaluate every 10 second of
    %simulated time) after half of pre-defined simulating time has passed
    if i>(T/dt)/2 && mod(i,10000)==0
        x1 =0+dx_mu/2:dx_mu:side-dx_mu/2; x2 =0+dx_mu/2:dx_mu:side-dx_mu/2;
        Inp=zeros(length(x1),length(x2));
        
        for idcell=1:N
            [X1,X2] = meshgrid(x1,x2);
            Inpu =MW(idcell)*PC_2DGaus([X1(:) X2(:)],mu(idcell,:),rate_noise);
            Inpu = reshape(Inpu,length(x2),length(x1));
            Inp(:,:)=Inp(:,:)+Inpu;
        end
        
    end
    
    if i>(T/dt)/2 && mod(i,10000)==0 &&  max(max(max(Inp)))<r_th
        no_plast_zones=1;
        fprintf('No more plasticity zones\n')
        break
        
    end
    
end

toc

% Plot weights and rat path, with a specific colormap for Fig. 5 and Fig. 6c
load('mwrcmap','mwrcmap')

MW_xy=reshape(MW(:),sqrt(N),sqrt(N));
pix_path=unique(round((x(1:i-7,:)+add_bord)/dx_mu),'rows');   % Select the pixels where the rat has been
for id_pix=1:length(pix_path)
    MW_xy(pix_path(id_pix,1),pix_path(id_pix,2))=0.01;
end
fi=figure;
imagesc(MW_xy')
set(gca,'YDir','normal')
xlim([ 1 sqrt(N)]);   set(gca,'XTick',[1 60 120 180]);     set(gca,'XTickLabel',[1 60 120 180])
ylim([ 1 sqrt(N)]);   set(gca,'YTick',[1 60 120 180]);     set(gca,'YTickLabel',[1 60 120 180])
set(fi,'Colormap',mwrcmap)
axis square


figure %not shown in the manuscript
x1 =0+dx_mu/2:dx_mu:side-dx_mu/2; x2 =0+dx_mu/2:dx_mu:side-dx_mu/2;

ExampInp=zeros(length(x1),length(x2));
for idpc=1:N
    [X1,X2] = meshgrid(x1,x2);
    Fplotw =(min(MW(idpc)+acetil(idpc),wmax))*PC_2DGaus([X1(:) X2(:)],mu(idpc,:),rate_noise);
    Fplotw = reshape(Fplotw,length(x2),length(x1));
    ExampInp=ExampInp+Fplotw;
end
imagesc(ExampInp)
set(gca,'YDir','normal')  ;    xlabel('m');   ylabel('m')
colorbar
axis square
title('Input(\Sigma w*r) to GC, available maze')


if use_rate
    
    % Generate grid cell firing map for a trial after training period - testing period:
    Tfam=(T/1)/dt;             % ms
    frfam=zeros(1,Tfam);
    tic
    for ifam=2:Tfam
        firrat_pc(:)=PC_2DGaus(x(ifam,:),mu,rate_noise);
        I(:)= MW'*firrat_pc(:);
        frfam(ifam)=frfam(ifam-1)+(dt/tau_r)*(-frfam(ifam-1)+SGtransf(I,0.12,0.5));
        
        if mod(ifam,round(Tfam/10))==0
            fprintf('Generating firing rate map. %d%% done. \n',round(100*ifam/Tfam))
        end
    end
    toc
    
    figure    %Fig. 6a ,Fig. 6d and Fig. 7d
    resoluc=400;
    x1 = linspace(0,side,resoluc*side); x2 = linspace(0,side,resoluc*side);
    RateSumMap=zeros(length(x1),length(x2));
    RateCountMap=zeros(length(x1),length(x2));
    tt=length(find(frfam(:)>0));
    BXS=ones(1,tt);
    BYS=ones(1,tt);
    for t=1:tt
        %a linha seguinte e similares carecem de simplificacao
        bx=max(1,resoluc*(round(resoluc*x(t,1)))/resoluc); %x coords for bin where the rat is at time t(the max is for no 0's as index)
        by=max(1,resoluc*(round(resoluc*x(t,2)))/resoluc); %y coords for bin where the rat is at time t(the max is for no 0's as index)
        RateSumMap(bx,by)=frfam(t)+RateSumMap(bx,by);
        RateCountMap(bx,by)=RateCountMap(bx,by)+1;
    end
    RateMapG=RateSumMap./RateCountMap;
    imagesc(RateMapG')
    set(gca,'YDir','normal')
    title('GC Rate map')
    axis square
        
    % GRIDSCORE COMPUTING
    
    GRIDMAPSTORE=zeros(sqrt(N_in),sqrt(N_in));
    
    resoluc=100;                  %For gridscore calculations we use resolution 100 (less time consuming)
    x1 = linspace(0,side,resoluc*side); x2 = linspace( 0,side,resoluc*side);
    RateSumMap=zeros(length(x1),length(x2));
    RateCountMap=zeros(length(x1),length(x2));
    tt=length(find(frfam(:)>0));
    for t=1:tt
        bx=max(1,resoluc*round(100*(x(t,1)))/100);
        by=max(1,resoluc*round(100*(x(t,2)))/100);
        RateSumMap(bx,by)=frfam(t)+RateSumMap(bx,by);
        RateCountMap(bx,by)=RateCountMap(bx,by)+1;
    end
    RateMap=RateSumMap./RateCountMap;
    RateMap(isnan(RateMap))=0;
    GRIDMAPSTORE(:,:)=RateMap;
end

score=0;

if use_rate
    GridMap=GRIDMAPSTORE(:,:);
    gaussWindow=fspecial('gaussian',7,10);    %Smooth the grid map to compute autocorrelogram
    SmoothedGridMap=conv2(GridMap,gaussWindow,'same');
    figure
    imagesc(SmoothedGridMap)
    set(gca,'YDir','normal')
    title('Grid Map using rate')
    axis square
    colormap;
    MaptoScore=SmoothedGridMap;
else
    GridMap=ExampInp;
    figure
    imagesc(GridMap )
    set(gca,'YDir','normal')
    title('Grid Map using synaptic input')
    axis square
    colormap;
    MaptoScore=GridMap;
    
    if x_haft
        red=fspecial('disk',90);
        red(red>0)=1;
        red=[red(1:89,1:89) red(1:89,91:181); red(91:181,1:89) red(91:181,91:181)];
        MaptoScore=MaptoScore.*red;
        figure;imagesc(MaptoScore);set(gca,'YDir','normal');axis square
    end
end
SpatialAutoCorr=xcorr2(MaptoScore);                       %Compute Autocorrelation Map
SpatialAutoCorrNorm=normxcorr2(MaptoScore,MaptoScore);    %Compute Normalized Autocorrelation Map

figure  % Fig. 6b
imagesc(SpatialAutoCorrNorm)
set(gca,'YDir','normal')
title('Normalized Autocorrelogram')
axis square

score=gridscore(SpatialAutoCorrNorm);

score
