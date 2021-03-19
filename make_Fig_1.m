%% make_Fig_1.m
%
% creates Figure 1 in the paper
%
% RW 18/3/21

clear all
close all

fresh_start = 0; % re-calculate data or load from a file?

if fresh_start

    % set constants
    alpha = 0.98; % stellar flux as fraction of present-day Earth value []
    nr    = 2^10;  % number of ensemble members / realizations []
    nt    = 5e3;  % number of timesteps for one realization []
    tau   = 2.5;  % linear relaxation timescale [My]
    sig_y = 0.7;  % distribution standard deviation []

    % initialize Climate and Solver objects
    climate = Climate();
    solver  = Solver(nr,nt,tau,sig_y);
    climate = climate.calculate_chi(alpha,solver);

    % solve the system
    solver = solver.solve_ensemble(climate);
    x_a    = climate.calculate_x(solver.y_a,alpha);

    % calculate numerical p.d.f.s
    [p_num, x_h] = solver.calculate_numerical_pdf(x_a);
    [q_num, y_h] = solver.calculate_numerical_pdf(solver.y_a);

    % convert numerical p.d.f.s to physical units
    fCO2_h = y_h*climate.fCO20; % CO2 molar concentration [ppmv]
    q_num  = q_num/climate.fCO20; % q(fCO2) [1/ppmv]
    T_h    = x_h + climate.T0; % temperature [K]

    % calculate analytical p.d.f.s
    chi  = climate.chi;
    dS   = climate.S0*(alpha-1);
    y    = @(x) exp((climate.a*x-dS)/climate.b);
    dydx = @(x) climate.a*y(x)/climate.b;
    %q    = @(y) solver.q(y,chi,sig_y);
    %p    = @(x) q(y(x)).*dydx(x);
    q_ana = solver.q(y_h,chi,sig_y);
    p_ana = solver.q(y(x_h),chi,sig_y).*dydx(x_h);

    % find a case where snowball transition occurred
    [it,ir_snow] = find(x_a<climate.xs & solver.y_a>2e-2,1);
    
    % save arrays for plotting
    t_a = solver.t_a;
    y_i = solver.y_a(:,ir_snow);
    x_i = x_a(:,ir_snow);
    ym_a = zeros(nt,1) + mean(solver.y_a(:,ir_snow))*climate.fCO20;
    xm_a = zeros(nt,1) + mean(x_a(:,ir_snow)) + climate.T0;
    
    % save some space
    clear solver x_a

else
    
    load fig_1_results.mat

end

% display fCO2 time series
subplot(2,2,1)
plot(t_a,y_i*climate.fCO20,'k'); hold on
plot(t_a,ym_a,'Color',[0.5 1 0.5]*0.75,'LineWidth',1.5)
xlabel('t [My]'); ylabel('f_{CO2} [ppmv]')
axis([0 500 0 2*climate.chi*climate.fCO20])

% display T time series
subplot(2,2,3)
plot(t_a,x_i+climate.T0,'k'); hold on
plot(t_a,xm_a,'Color',[0.5 1 0.5]*0.85,'LineWidth',1.5)
xlabel('t [My]'); ylabel('T [K]')
h = patch([0 500 500 0], [278 278 280 280],'c','LineStyle','none');
set(h,'facealpha',.2)
text(25,279,'Snowball transition region')
axis([0 500 278 290])

% display fCO2 pdf
subplot(2,2,2)
plot(q_num,fCO2_h,'k',q_ana/climate.fCO20,fCO2_h,'r--')
xlabel('q(f_{CO2}) [1/ppmv]'); ylabel('f_{CO2} [ppmv]')
axis([0 3e-3 0 2*climate.chi*climate.fCO20])

% display T pdf
subplot(2,2,4)
semilogx(p_num,T_h,'k',p_ana,T_h,'r--')
xlabel('p(T) [1/K]'); ylabel('T [K]')
axis([1e-5 1 278 292])
xticks([1e-5 1e-4 1e-3 1e-2 1e-1 1])
