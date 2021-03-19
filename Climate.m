%% Climate.m
%
% Contains general climate constants and routines
%
% RW 18/3/21

classdef Climate
    
    properties (SetAccess = public)
        
        F0 % Earth present-day solar constant [W/m2]
        A0 % Earth present-day albedo []
        S0 % Earth present-day absorbed solar radiation [W/m2]
        T0 % Earth pre-industrial surface temperature [K]
        fCO20 % Earth pre-industrial CO2 molar concentration [ppmv]
        a % OLR linearization parameter [W/m2/K]
        b % OLR linearization parameter [W/m2]
        xs % Snowball transition temperature deviation from 288 K [K]
        chi % y probability distribution parameter [mol/mol]
        log_ym % mean of log(y)
        
    end
    
    methods
        
        function self = Climate()
            
            % scalars
            self.F0 = 1366;
            self.A0 = 0.3;
            self.S0 = self.F0*(1-self.A0)/4;
            self.T0 = 288;
            self.fCO20 = 280;
            self.a = 2;
            self.b = 5.35;
            self.xs = -8;

        end
        
        function self = calculate_chi(self,alpha,solver)
            % calculate chi for a given stellar flux

            % calculate <log y> assuming <x>=0
            self.log_ym = self.S0*(1 - alpha)/self.b;
            
            % first guess for chi
            chi_ini = exp(self.log_ym);
            
            % calculate actual chi value
            if(chi_ini/solver.sig_y>15)
                self.chi = chi_ini;
            else
                Phi = self.calculate_Phi(solver,self.log_ym);
                self.chi = fzero(Phi,chi_ini);
                self.chi = max(0,self.chi);
            end
            
        end
        
        function Phi = calculate_Phi(~,solver,log_ym)
            % calculate Phi(chi) anon. fn. (eqn. 12 in text)

            % set up integration grid to always capture p.d.f. with sufficient
            % resolution
            ya_min = max(1e-8,exp(log_ym) - solver.sig_y*15);
            ya_max = exp(log_ym) + solver.sig_y*15;
            ya = linspace(ya_min,ya_max,1e4);
            Phi = @(chi)  log_ym - trapz(ya,log(ya).*solver.q(ya,chi,solver.sig_y));
        end
        
        function x = calculate_x(self,y,alpha)
            % calculate x (warm branch) for a given y

            deltaS = -self.S0*(1 - alpha);
            x = (deltaS + self.b*log(y))/self.a;

        end
        
        function alpha = calculate_alpha(~,t)
            % calculate alpha for a given t
        
            alpha = 1./(1 + (2/5)*(1 - t/4.5e3));

        end
        
    end
end

