%% make_Fig_2.m
%
% creates Figure 2 in the paper
%
% RW 18/3/21

clear all
close all

fresh_start =0; % re-calculate data or load from a file?

if fresh_start
    
    % set constants
    nr  = 1; % number of realizations []
    nt  = 5.0e4; % number of timesteps for one realization []
    ns  = 4; % number of sigma values to use []
    tau = 2.5; % linear relaxation timescale [My]
    
    % initialize Climate object
    climate = Climate();

    % initialize x, y, sigma and color scheme arrays
    sig_y_a = [0.1 1 10 100];
    cvec    = [0.1 0.1 0.0; 0.3 0.1 0.0; 0.6 0.1 0.0; 0.9 0.1 0.0];
    x_a     = zeros(ns,nt);
    y_a     = zeros(ns,nt);
    
    for is=1:ns
        
        % assign std. dev. value
        sig_y = sig_y_a(is);
        
        % initialize Solver object
        solver = Solver(nr,nt,tau,sig_y);

        % calculate alpha vs. time
        alpha_a = climate.calculate_alpha(solver.t_a);
        
        % calculate initial chi value
        climate = climate.calculate_chi(alpha_a(1),solver);
        
        % initial conditions
        y = climate.chi; % a good starting approx. when chi >> sigma_y 
        x = 0;
        
        % start time iteration
        for it=1:nt
            
            % update chi value
            climate = climate.calculate_chi(alpha_a(it),solver);

            % update y
            y = solver.single_step(y,climate);
            
            % stop run if temperature drops below Snowball threshold
            if(x<climate.xs)
                x_mean(is) = mean(x_a(is,1:it));
                x_a(is,it:end) = NaN;
                break
            end
            
            % update x
            x = climate.calculate_x(y,alpha_a(it));
            
            % save x and y to 2-d arrays
            x_a(is,it) = x;
            y_a(is,it) = y;
            
        end
        
        is
        
    end
    
    % calculate time as Gy before present
    t_Gya = 4.5 - solver.t_a/1e3;
            
else
    
    load fig_2_results.mat
    
end

% used to create 2nd x-axis for plot (this part done in Adobe
% Illustrator)
alpha_plt = climate.calculate_alpha([1.5 2 2.5 3 3.5 4 4.5 5]*1e3)

for is=1:4
    
    % display fCO2 vs. time
    h1 = subplot(2,1,1)
    semilogy(t_Gya,y_a(is,:)*climate.fCO20,'Color',cvec(is,:)); hold on
    xlabel('time [Gy before present]')
    ylabel('f_{CO2} [ppmv]')
    axis([-0.5 3 1 1e6])
    set(h1, 'Xdir', 'reverse')
    
    % display  T vs. time
    h2 = subplot(2,1,2)
    plot(t_Gya,x_a(is,:)+climate.T0,'Color',cvec(is,:)); hold on
    xlabel('time [Gy before present]')
    ylabel('T [K]')
    axis([-0.5 3 278 295])
    set(h2, 'Xdir', 'reverse')
    
end

% add legend and Snowball transition region shading

subplot(2,1,1)
legend('\sigma_y = 0.1','\sigma_y = 1','\sigma_y = 10','\sigma_y = 100')

subplot(2,1,2)
h = patch([-0.5 3 3 -0.5], [275 275 280 280],'c','LineStyle','none');
set(h,'facealpha',.15)
text(2.9,277.5,'Snowball transition region')
axis([-0.5 3 275 295])
