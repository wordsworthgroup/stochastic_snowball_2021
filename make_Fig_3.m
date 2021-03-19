%% make_Fig_3.m
%
% creates Figure 3 in the paper
%
% RW 18/3/21

clear all
close all

fresh_start = 0; % re-calculate data or load from a file?

if fresh_start
    
    nr  = 1; % number of realizations []
    nt  = 5.0e4; % number of timesteps for one realization []
    tau = 2.5; % linear relaxation timescale [My]
    
    % initialize Climate object
    climate = Climate();
    
    % set up sig_y_a and other arrays
    sig_y_a    = logspace(-1,2,2^5);
    alpha_trans = zeros(length(sig_y_a),1);
    t_trans     = zeros(length(sig_y_a),1);
    chi_a       = zeros(1,nt);
    
    % loop over std. dev. values
    for is=1:length(sig_y_a)
        
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

        for it=1:nt
            
            % update chi (does not vary with realization)
            climate = climate.calculate_chi(alpha_a(it),solver);
            
            % get chi value (does not vary with realization)
            if(is==1)
                chi_a(it) = climate.chi;
            end
            
            % update y
            y = solver.single_step(y,climate);
            
            % stop run if temperature drops below Snowball threshold
            if(x<climate.xs)
                alpha_trans(is) = alpha_a(it);
                t_trans(is)     = solver.t_a(it);
                break
            end
            
            % update x
            x = climate.calculate_x(y,alpha_a(it));
            
        end
        
        is
        
    end
    
    % calculate sigma_f
    sig_f = sig_y_a*climate.fCO20;
    
else
    load fig_3_results.mat
end

% display Snowball transition time vs. sigma_y 
% also plot chi vs. time
h1 = subplot(1,1,1);
semilogx(chi_a*climate.fCO20,4.5-solver.t_a/1e3,'r'); hold on
plot(chi_a*climate.fCO20/2,4.5-solver.t_a/1e3,'r--')
plot(chi_a*climate.fCO20/4,4.5-solver.t_a/1e3,'r:')
scatter(sig_f,4.5-t_trans/1e3,50,'k.');
grid on
xlabel('\sigma_y [ppmv]')
ylabel('time of Snowball transition [Gy before present]')
set(h1, 'Ydir', 'reverse')
axis([1e1 1e5 -.5 2])
legend('\chi','\chi/2','\chi/4','model output')
