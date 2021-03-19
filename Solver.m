%% Solver.m
%
% Contains core routines to solve the stochastic evolution equation
% RW 18/3/21

classdef Solver
    
    properties (SetAccess = public)
        
        nr % number of realizations []
        nt % number of timesteps for one realization []
        tau % linear relaxation timescale [My]
        dt % timestep [My]
        t_a % time array
        y_a % y array
        sig_y % std. dev. of y []
        g % noise ampitude [1/My^0.5]
        Q % normalization constant []
        q % probability distribution for y []
        
    end
    
    methods
        
        function self = Solver(nr,nt,tau,sig_y)
            
            % scalars
            self.dt = 0.1;
            self.nr = nr;
            self.nt = nt;
            self.tau = tau;
            self.sig_y = sig_y;
            self.g = sqrt(2/self.tau)*self.sig_y;

            % arrays
            self.y_a = zeros(self.nt,self.nr); % y array []
            self.t_a = (1:self.nt)*self.dt; % time array [My]
            
            % Q-anon functions [eqns. 6 and 7 in text]
            self.Q = @(chi,sig_y) (2/(sig_y*sqrt(2*pi)))/(1 - erf(-chi/(sqrt(2)*sig_y)));
            self.q = @(y,chi,sig_y) self.Q(chi,sig_y)*exp(-(y-chi).^2/(2*sig_y^2));                
            
        end
        
        function y = single_step(self,y,climate)
            % do a single stochastic model step

            % calculate Wiener step and update y
            dW = sqrt(self.dt)*randn(1,1);
            dy = -(y-climate.chi)*self.dt/self.tau + self.g*dW;
            y = y + dy;
            
            % boundary condition to avoid unphysical -ve y
            if(y<0)
                y=-y;
            end
            
        end
        
        function self = solve_ensemble(self,climate)
            % solve for an ensemble of simulations
            for ir=1:self.nr
                
                % initialize y 
                % (gives statistically steady state for chi>>sig_y only)
                y = climate.chi + self.sig_y*randn(1,1);
                
                % boundary condition to avoid unphysical -ve y
                if(y<0)
                    y=-y;
                end
                
                % solve using Euler method
                for it=1:self.nt
                    y = self.single_step(y,climate);
                    self.y_a(it,ir) = y;
                end
                
            end
        end
        
        function [f_num, x_bins] = calculate_numerical_pdf(~,x_a)
            % calculate numerical pdf
        
            hist_tmp = histogram(x_a,150);
            f_num    = hist_tmp.Values;
            x_bins   = hist_tmp.BinEdges(1:end-1) + hist_tmp.BinWidth/2;
            f_num    = f_num/trapz(x_bins,f_num);

        end
        
    end
end
