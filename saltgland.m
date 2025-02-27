%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function saltgland.m
% Melissa H Mai
%
% Solves the secreting salt gland model using Newton's method
%
% INPUTS (mandatory)
% rp        Fracture sizes (m). Can be passed as a vector.
% Jgc       Maximum salt transporter rate (kg/m^2/s).
%
%
% INPUTS (optional)
% Optional inputs should be passed as 'Name', Value pairs when calling the
% function (ie, saltgland(rf, Jgc, 'cstar', 0.8, 'r', 10e-6). Input names
% are not case-sensitive.
% 
% 'J', 'JFunc'          String describing the transporter activity profile:
%                           'sigmoid': Jgc./(1 + exp(-h*(cstar-c+cg)))
%                           'exponential': Jgc.*exp(-2*log(10)*(c-cg)/cstar)
%                           'linear': max(Jgc(1)*(1-(c-cg)/cstar),0)
%                           Default: 'sigmoid'
% 'cstar', 'ccrit'      Threshold concentration gradient for transporter
%                           activity (mM). 
%                           Default: 1400
% 'cg'                  Concentration within the gland (mM).
%                           Default: 100
% 'h'                   Profile steepness for sigmoid activity profiles.
%                           Default: 100
% 'lp'                  Cuticle thickness (m)
%                           Default: 1e-6
% 'D'                   Diffusion coefficient of sodium (m^2/s)
%                           Default: 1.6e-9
% 'r'                   Gland radius (m)
%                           Default: 1.5e-5
% 'eta'                 Viscosity of water (kg/m/s)
%                           Default: 1e-3
% 'n'                   Number of cracks
%                           Default: 1
% 'Lmp'                 Membrane permeability (m^2 s/kg)
%                           Default: 1e-12
% 'cb'                  Brine concentration (mM). If Pb is defined, 
%                           cb = (Pb - psi_b)/RT
%                           Default: 2000
% 'psi_g'               Gland water potential (kg/m/s^2, Pa)
%                           Default: -1.5e6
% 'psi_b'               Brine water potential (kg/m/s^2, Pa). If cb is
%                           defined, psi_b = Pb - RT*cb
%                           Default: -5e6
% 'T', 'RT'             RT (thermal energy scale) (kg m^2/s^2)
%                           Default: 298*8.314
%
%
% OUTPUT
% The output is a structure `output` with the following fields:
% c                     Concentration within the chamber (dimensionless,
%                           scaled by cb)
% cb                    Brine concentration (mM)
% Pc                    Chamber hydrostatic pressure (MPa)
% psi_c                 Chamber water potential (MPa)
% Jtot                  Net water loss (mol/s)
% Deff                  Calculated effective diffusion across the cuticle (m^3/s)
% Lcp                   Calculated conductance across the cuticle (m^4 s/kg)
% JD                    Diffusive component of salt flux (mol/s)
% JP                    Advective component of salt flux (mol/s)
% Jsalt                 Net salt efflux (mol/s)
% Jfunc                 Function handle for the activity profile
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function output = saltgland(rf, Jgc, varargin)
    % Process inputs
    for i = 1:2:length(varargin)
        this = varargin{i+1};
        switch lower(varargin{i})
            case {'cstar','ccrit'}
                cstar = this;
            case {'jfunc', 'j'}
                Jfunc = this;
            case 'lp'
                lp = this;
            case 'd'
                D = this;
            case 'r'
                r = this;
            case 'eta'
                eta = this;
            case 'n'
                n = this;
            case 'lmp'
                lmp = this;
            case 'cb'
                cb = this;
            case 'psi_g'
                psi_g = this;
            case 'psi_b'
                psi_b = this;
            case {'t', 'rt'}
                RT = this;
            case 'h'
                h = this;
            case 'cg'
                cg = this;
        end
    end

    Pb = 1.01325e5;

    % Default values
    if ~exist('cstar','var'); cstar = 1400; end % mM
    if ~exist('cg', 'var'); cg = 100; end % mM (mmol/kg)
    if ~exist('Jfunc','var'); Jfunc = 'sigmoid'; end
    if ~exist('lp', 'var'); lp = 1e-6; end % m
    if ~exist('D', 'var'); D = 1.6e-9; end % m^2/s
    if ~exist('r', 'var'); r = 1.5e-5; end % m
    if ~exist('eta', 'var'); eta = 1e-3; end % kg/m/s
    if ~exist('n', 'var'); n = 3; end % Dimensionless
    if ~exist('lmp', 'var'); lmp = 1e-12; end % m^2 s / kg
    if ~exist('RT', 'var'); RT = 298*8.314; end % kg m^2 / s^2
    if ~exist('psi_g', 'var'); psi_g = -1.5e6; end % kg / m / s^2
    if ~exist('cb','var'); cb = 2000; end % mmol / kg = mol / m^3
    if ~exist('psi_b', 'var') % kg / m / s^2
        if exist('cb','var')
            psi_b = Pb - RT*cb;
        else
            psi_b = -5e6;
        end
    end
    if ~exist('h', 'var'); h = 100; end

    if cstar > 1
        cstar = cstar/cb;
    end

    % Brine concentration
    p0 = Pb;
    cb = (Pb - psi_b)/RT;

    % Gland concentration
    cg = cg/ cb;
    
    % Dimensional analysis
    eta = eta*D/(p0*lp^2);
    r = r/lp;
    lmp = lmp*p0*lp/D;
    psi_g = psi_g/p0;
    psi_b = psi_b/p0;
    Pb = Pb/p0;
    RT = RT*cb/p0;
    
    rf = rf/lp;
    Jgc = Jgc*lp*pi*r^2./(D*cb); % Also made extensive
    
    % Effective diffusion
    Deff = 4*rf + n.*pi*rf.^2;

    % Cuticle resistance
    R = (eta./n).*(8./(pi*rf.^4) + 3./(rf.^3));

    % Water molar volume
    vbar = 1.807e-5 * cb;

    % Permeability of the cuticle
    Lcp = 1./R;
    
    % Permeability of the membrane
    Lmp = lmp*pi*r.^2;
    
    % Scaled
    xi = Lmp./Lcp;

    %%
    if strcmp(Jfunc,'sigmoid')
        J = @(c) max(Jgc./(1 + exp(-h*(cstar-c+cg)))-0.01*Jgc,0);
    elseif strcmp(Jfunc,'linear')
        J = @(c) max(Jgc(1)*(1-(c-cg)/cstar),0);
    elseif strcmp(Jfunc,'exponential')
        J = @(c) max(Jgc.*exp(-2*log(10)*(c-cg)/cstar)-0.01*Jgc,0);
    end

    % First order derivative of transporter function
    Jp = @(c) (J(c+.01)-J(c-.01)/.02);

    % Coefficients for psi_c = alpha*cc + beta
    alpha = - RT./(1 + xi);
    beta = (xi.*psi_g + Pb)./(1+xi);

    % Equation to solve for Newton's method
    F = @(c)(J(c)+alpha*Lmp.*c.^2-Deff.*c+Lmp.*(beta-psi_g)+Deff);
    Fp = @(c) (Jp(c) + 2.*alpha.*Lmp.*c - Deff);

    % Initialize
    cc = zeros(size(Jgc.*xi));
    k = 1;
    niter = 200;
    epsilon = 1e-9;

    % Newton's method
    while (sum(sum(abs(F(cc))))>epsilon && k < niter)
        cc = cc - F(cc)./Fp(cc);
        k = k+1;
    end

    %% Finish up
    psi_c = alpha.*cc + beta;

    % Final salt flux
    Jsalt = J(cc);

    % Remove vanishingly small numbers
    Jsalt(Jsalt<1e-8) = 0;

    % Diffusive Flux
    JD = -Deff.*(1-cc);
    
    % Cavity water pressure
    Pc = psi_c + cc.*RT;
    
    % Pressure Flux
    JP = -Lcp.*cc.*(Pb - Pc);

    % Total water flow, mol/s
    Jtot = Lmp*(psi_g - psi_c)/vbar;
    
    % Write to output
    output = struct('c', cc, 'Pc', Pc*p0*1e-6, 'psi_c', psi_c*p0*1e-6, ...
        'Jtot', Jtot*D*lp*cb, 'Deff', Deff*D*lp, 'JD', JD*D*lp*cb, ...
        'Jsalt', Jsalt*D*lp*cb, 'Jfunc', J,...
        'JP', JP*D*lp*cb, 'cb', cb, 'Lcp', Lcp*D*lp/p0);
    
end
