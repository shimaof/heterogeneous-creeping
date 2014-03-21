function solver_heterogeneous_creeping(model,test)

if nargin<1 model =1; test = 2; end  % by default, creeping model in the creeping simulation


clc;  % clear the command window
close all % close all figures

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Examples and comparisons with the n-populations model.
% For the input parameters
%     model = 1: the creeping model
%     model = 2: the n-populations model
%     test = 1: the overtaking simulation
%     test = 2: the creeping simulation
% For other parameters in the code
%     rho1, rho2 are densities of first and second vehicle class
%     rm1, rm2 are effective jam densities of rho1 and rho2
%     vm1, vm2 are the maximum traffic veloicty of rho1 and rho2
%------------------------------------------------
% This solver is based on the finite volume Godunov method, and have
% applied the HLL approximate Riemann solver to find the numerical class.
% In this simulation, first vehicle class rho1 is assumed to be small
% vehicles, while the second vehicle class is larger.

%--------------------------------------------------
% Aug 29 2013
% Shimao Fan
% University of Illinois at Urbana-Champaign
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

% parameter of effective jam density and maximim velocity
switch model
    case 1
       model_name = 'Creeping model';
       rm1 = 1.5; rm2 = 1.0;
%        vm1 = 1.0; vm2 = 2.0;
       vm1 = 1.5; vm2 = 1.5;
    case 2
       model_name = 'n-populations model';
       rm1 = 1.5; rm2 = 1.5;
%        rm1 = 1.0; rm2 = 1.0;
       vm1 = 1.5; vm2 = 1.0;
end
% ----------------------------------
% class 1 is small vehicle
% class 2 is large vehicle
%----------------------------------
rm = max(rm1,rm2);

tfinal = 250;
% tfinal = 50;
len = 50;
x = linspace(0,len,1000);   % space step
dx = x(2)-x(1);
% time step
k = (x(2)-x(1))/vm1;

lambda = k/dx;
dt = k;
t = 0:dt:tfinal;
M = length(t);   % length of time vector
N = length(x);   % length of space vector

%==========================================
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% boundary: traffic light, red
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
switch test
    case 1  % overtaking

     U(1,1:ceil(1*N/50)) = 0.00000000;
     U(1,ceil(1*N/50)+1:ceil(8*N/50)) = 0.8;
     U(1,ceil(8*N/50)+1:N) = 0;

     U(2,1:ceil(8*N/50)) = eps;
     U(2,ceil(8*N/50)+1:ceil(15*N/50)) = 0.8;
     U(2,ceil(15*N/50)+1:N) = eps;
     
     d1l = U(1,1)+0*t;
     d2l = U(2,1)+0*t;

     d1r = U(1,end)+0*t;
     d2r = U(2,end)+0*t;

     
    case 2 % creeping
      U(1,1:ceil(1*N/50)) = 0;
      U(1,ceil(1*N/50)+1:ceil(15*N/50)) = .7;
      U(1,ceil(15*N/50)+1:N) = 0;

      U(2,1:ceil(15*N/50)) = eps;
      U(2,ceil(15*N/50)+1:N) = 0.7;
     
      d1l = U(1,1)+0*t;
      d2l = U(2,1)+0*t;

      d1r = .89999+0*t;
      d2r = 0.6+0*t;
%       d1r = .5+0*t;
%       d2r = 1.0+0*t;

end
%----------------------------------------
U_0 = U;    % initial condition
inter = 250;   % for printing every 250 steps
% inter = 50;



%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%        The main algorithm
%/////////////////////////////////////

for n = 0:M-1
    % the boundary condition
    lbc = [d1l(n+1);d2l(n+1)];
    rbc = [d1r(n+1);d2r(n+1)];
    % insert the boundary condition
    
    Up1=[U(:,2:end),rbc];    % right cell bounds
    Um1=[lbc,U(:,1:end-1)];  % left cell bounds

    
    
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    % phase transition happens
    %////////////////////////////
    if model ==1
      rl = Um1(1,:) + Um1(2,:);   % total density on left
      r = U(1,:) + U(2,:);        % total density in middle
      rr = Up1(1,:) + Up1(2,:);   % total density on right 

      case1 = rl<=rm2;        % left state in the domain D1
      case2 = r<=rm2;         % middle state in the domain D1
      case3 = rr<=rm2;        % right state in the domain D1

      id_l = find(case1-case2); % find all phase changes ul and u
      id_r = find(case2-case3); % find all phase changes u and ur
    
      idl = case1(id_l)-case2(id_l);% -1 means left large density
      idr = case2(id_r)-case3(id_r);% +1 mean left state small
    
% 
      for i = 1:length(idr)-1
        if idr(i)<0
           U(1,id_r(i)) = rm2-U(2,id_r(i));   % insert a intermediate state
        else
           Up1(1,id_r(i)) = rm2-Up1(2,id_r(i));   % insert a intermediate state
        end
      end
    end
    %============================================
    % updating with godunov methods/approximate riemann solver
    U = U+lambda*(HLL(Um1,U,vm1,vm2,rm1,rm2)-HLL(U,Up1,vm1,vm2,rm1,rm2));
    %============================================
    
    
    % \\\\\\\\\\\\\\\
    % plotting part
    %/////////////////
    if mod(n,inter) ==0
    subplot(2,1,1),
    plot(x(2:1:end-5),U_0(1,2:1:end-5),'--','color',[0.8,0,0],'linewidth',2), hold on
    plot(x(2:end-5),U(1,2:end-5),'-','color',[.6,.6,.6],'linewidth',6)
%     title(sprintf('%s, t=%3.0f',model_name,n*dt),'fontsize',25)
    title(sprintf('%s',model_name),'fontsize',25)
    axis([0 len 0 1.05*rm])
    h = legend('Initial state','Density of small vehicles');
    set(h,'Location','NorthWest')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',25)    
    set(gca,'xtick',[])     
    set(gca,'position',[.10 .53 .885 .41])
    ylabel('First vehicle class \rho_1')
    
    %****************
    subplot(2,1,2),
    %****************
    plot(x(2:1:end-5),U_0(2,2:1:end-5),'--','color',[0.8,0,0],'linewidth',2), hold on
    plot(x(2:end-5),U(2,2:end-5),'-.','color',[0,0,.1],'linewidth',6)
    h = legend('Initial state','Density of large vehicles');
    set(h,'Location','NorthWest')
    set(gca,'linewidth',1.5)
    set(gca,'fontsize',25)
    axis([0 len 0 1.05*rm])
    ylabel('Second vehicle class \rho_2')
    xlabel('x')
   
    
    set(gca,'position',[.10 .1 .885 .41])
    res = 800;
    set(gcf,'paperpositionmode','auto')
    set(gcf,'position',[10  50 res res*.8])
     
%      print figure
%      filename_save = sprintf('fig_multi_class_%1.0f_%s_%1.0f',test,model_name(1:3),ceil(n/inter));
%      print(gcf,'-dpng',filename_save,'-r290') 
%      print(gcf,'-depsc',filename_save,'-r290','-painters')
%      fixPSlinestyle([filename_save,'.eps'])
    end
    drawnow
end





%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% subfunctionals
%///////////////////////////////////////////////////////


function F = HLL(Ul,Ur,vm1,vm2,rm1,rm2) % the HLL Riemann solver
    % where v contains boundary condition of velocity
    % eigenvalues:
    eigns1 = lambda(Ul,vm1,vm2,rm1,rm2);  % left bound
    eigns2 = lambda(Ur,vm1,vm2,rm1,rm2);  % right bound
    % another way to define c1 and c2
    c1 = min(eigns1(1,:),eigns2(1,:));
    c2 = max(eigns1(2,:),eigns2(2,:));
    case1 = c1>=0;
    case2 = c1<0 & c2>=0;
    case3 = c2<0;
    % left flux
    Fl = max(flux(Ul,vm1,vm2,rm1,rm2),0);  % plug into left boundary velocity information
    % right flux
    Fr = max(flux(Ur,vm1,vm2,rm1,rm2),0);
    %-----------------------------------------------
    % intermediate flux
    F_hll(1,:) = (c2.*Fl(1,:)-c1.*Fr(1,:)+(Ur(1,:)-Ul(1,:)).*c1.*c2)./...
            (c2-c1);
    F_hll(2,:) = (c2.*Fl(2,:)-c1.*Fr(2,:)+(Ur(2,:)-Ul(2,:)).*c1.*c2)./...
            (c2-c1);
    F(1,:) = case1.*Fl(1,:)+case2.*F_hll(1,:)+case3.*Fr(1,:);
    F(2,:) = case1.*Fl(2,:)+case2.*F_hll(2,:)+case3.*Fr(2,:);

%--------------------------------------------------------------------------   
function y = lambda(U,vm1,vm2,rm1,rm2)  % U = (density, w), calculate eigenvalues
       r = U(1,:)+U(2,:);    % sum
       u1 = vel(r,vm1,rm1);
       u2 = max(vel(r,vm2,rm2),0);
       dv1 = DiffVel(r,vm1,rm1);
       indx = find(r>=rm2);
       dv2 = DiffVel(r,vm2,rm2);
       dv2(indx) = 0;
       % characteristic speed
       c1 = U(1,:).*dv1;
       c2 = U(2,:).*dv2;
       q1 = c1+u1; q2 = u2+c2;
       delta = sqrt((q1-q2).^2+4*c1.*c2);
       
       y(1,:) = .5*(q1+q2-delta);   % slower chareateritic field
       y(2,:) = .5*(q1+q2+delta);   % faster characteristic field
%-----------------------------------------------
function y = flux(U,vm1,vm2,rm1,rm2)  % define flux function        
      u1 = vel(U(1,:)+U(2,:),vm1,rm1);
      u2 = vel(U(1,:)+U(2,:),vm2,rm2);
      y(1,:) = U(1,:).*u1;
      y(2,:) = U(2,:).*u2;
%-----------------------------------------------
function y = DiffVel(r,vm,rm)% define derivative of velocity        
      y = (vel(r+eps,vm,rm)-vel(r,vm,rm))/eps;
%-------------------------------------------------------
function y = vel(rho,vm,rhom)  % define velocity functions 
     y = vm*(1-(rho/rhom));   % which is the linear form 
         
%\\\\\\\\\\\\\\\\\\\\\\\\THE END OF CODE\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
