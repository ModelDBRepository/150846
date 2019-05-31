% This code computes the mechanism for Plasticity Rule demonstration 
% with canonical neural operations
% Luisa Castro, FCUP
% luisa.castro@fc.up.pt

clear all
close all

T=10000;            % [ms] simulation duration =10s
dt=1;               % [ms]time step for each iteration

N=3;                % number of intermediate neurons in the circuit

r=zeros(N,T);       % storage of firing rate values
tau_r=100;          % [ms] time constant for neurons firing rate

dw_e=zeros(1,T);    % storage weight modification amounts for excitatory synapse, from N1 to GC 
dw_i=zeros(1,T);    % storage weight modification amounts for inhibitory synapse, from N2 to GC

alpha=0.01;         % scaling parameter for synapse modification
a=0.02;             % slope parameter for sigmoid transfer function
b1=1/4;             % inflexion points for sigmoid transfer function
b3=2/4;
b2=3/4;

rPC=linspace(0,1,T);   % input rates taking every value in the normalized firing rate range, smoothly

for i=2:T
    
    r(1,i) = r(1,i-1)  +  (dt/tau_r) * (  -r(1,i-1)  +  SGtransf(rPC(i)        ,a,b1  ) );
    r(2,i) = r(2,i-1)  +  (dt/tau_r) * (  -r(2,i-1)  +  SGtransf(rPC(i)       ,a,b2  ) );
    r(3,i) = r(3,i-1)  +  (dt/tau_r) * (  -r(3,i-1)  +  SGtransf(rPC(i)-r(2,i),a,b3  ) );
    
    dw_e(i) = alpha*r(1,i);
    dw_i(i) = alpha*r(3,i);
    
end

figure   				% Fig 4
subplot(411)       
plot(rPC,r(1,:),'LineWidth', 2,'color','b');    
set(gca,'XTick',[0 b1 b3 b2 1]);    set(gca,'XTickLabel',{'','','','',''})
ylabel('Rate');                
ylim([-0.05 1.05]);                 set(gca,'YTick',[0 1]); set(gca,'YTickLabel',{'0','1'}); 
hold on
plot(rPC,r(3,:),'LineWidth', 2,'color','g');
plot(rPC,r(2,:),'LineWidth', 2,'color','r'); 
h=legend('E','I_1','I_2',3);        set(h,'box','off','Location','NorthWest')
box off

subplot(412)
plot(rPC,dw_e,'r','LineWidth', 2,'color','b'); 
set(gca,'YTick',[0 0.01]);          set(gca,'XTick',0:0.5:1); set(gca,'XTickLabel',{'','',''})
ylabel('\Deltaw_+');       
ylim([-0.0005 0.0105]);             set(gca,'YTickLabel',{'',''})
box off
       
subplot(413)
plot(rPC,dw_i,'b','LineWidth', 2,'color','g'); 
set(gca,'YTick',[0 0.01]);          set(gca,'XTick',0:0.5:1); set(gca,'XTickLabel',{'','',''})
ylabel('\Deltaw_-');     
ylim([-0.0005 0.0105]);             set(gca,'YTickLabel',{'',''})
box off

subplot(414)
plot(rPC,(dw_e-dw_i),'LineWidth', 2,'color','b'); 
xlabel('PC rate');          
xlim([0 1]);           
set(gca,'XTick',[0 b1 b3 b2 1]);    set(gca,'XTickLabel',{'0','b1','b3','b2','1'})
ylabel ('\Deltaw ');     
ylim([-0.0005 0.0105]);             set(gca,'YTick',[0 0.01]);  set(gca,'YTickLabel',{'',''})
box off