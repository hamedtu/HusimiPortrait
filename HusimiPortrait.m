clear; close all;
%%
A = load("TS_eigenvectors_0.92_0.40_1.90_10000_64_real.dat");
B = load("TS_eigenvectors_0.92_0.40_1.90_10000_64_imag.dat");

FloquetMatrix = A + 1i * B;
N = 10000;
alpha = 0.92;
mu = 0.4;
w = 1.9;

[eta, I] = Coherence(FloquetMatrix, N);

% Figure 3
index = [1, 110, 194, 276, 768, 972, 1415, 1673];
Husimiplot(FloquetMatrix, I(index), [alpha mu w], N, 250);
% Figure 6
index2 = 5705;
Husimiplot(FloquetMatrix, I(index2), [alpha mu w], N, 250);
%% Function Definitions
%%
% Function: pSpacePortrait
%
% Description:
%   Generates and plots a Poincare surface of section (phase-space portrait)
%   by integrating a two-dimensional (4-dimensional) dynamical system.
%
% Inputs:
%   - param1: Meanfield parameter alpha
%   - param2: Driving amplitude
%   - param3: Driving frequency
%
function pSpacePortrait(param1, param2, param3)
    % Define grid for p and phi
    p = linspace(-0.9, 0.9, 35);
    phi = linspace(-3.1, 3.1, 35);
    
    % Define the period of the driving
    T = 2*pi/param3;
    
    % Define the equations of motion
    dgl = @(t, x, param1, param2, param3) [ -sqrt(1 - x(1)^2) * sin(x(2));
                                           2.0*param1*x(1) + x(1)/sqrt(1 - x(1)^2) * cos(x(2)) + 2*param2*cos(param3*t)];
    
    % Define time limits and tolerance for integration
    tLim = 0:T:200*T;
    options = odeset('RelTol', 1e-7, 'AbsTol', 1e-10);
    
    hold on
    % Loop over the grid and compute trajectories
    for k = 1:length(p)
        for l = 1:length(phi)
            [t, y] = ode45(@(t, y) dgl(t, y, param1, param2, param3), tLim, [p(k); phi(l)], options);
            % Adjust phase angle
            Phi = y(:, 2) / pi;
            Phi = Phi + 1;
            Phi = mod(Phi, 2);
            Phi = Phi - 1;
            % Plot the trajectory as a scatter plot
            scatter(Phi, y(:, 1), 1.0, [0, 0, 0], 'filled');
        end
    end
    hold off
    xlabel('$\varphi / \pi$', 'Interpreter', 'latex')
    ylabel('$p$', 'Interpreter', 'latex')
end

%%
% Function: Husimiplot
%
% Description:
%   Computes and plots the Husimi distributions for one or more Floquet states.
%   Each state's Husimi matrix is normalized by its maximum value before
%   summing (for visibility when considering several states). The function uses a specified grid in phase space (p and phi) and
%   overlays the phase space portrait.
%
% Inputs:
%   - V: Matrix containing Floquet states (each column is a state).
%   - V_ind: Vector of indices referring to the states to be processed.
%   - param: A three-element vector [param1, param2, param3] where:
%            param1 is the mean-field parameter, param2 the driving amplitude,
%            and param3 the driving frequency.
%   - N: Particle number.
%   - gridSize: Resolution of the phase space grid.
%
function Husimiplot(V, V_ind, param, N, gridSize)
    % Initialize the phase space grid
    p = linspace(0, 1, gridSize);
    phi = linspace(-pi, pi, gridSize);
    Hussimi = zeros(gridSize, gridSize);
    
    % Loop over each index in V_ind to compute each Husimi matrix
    for j = 1:length(V_ind)
        v = V(:, V_ind(j));
        HussimiMatrix = zeros(gridSize, gridSize);
        
        % Compute the Husimi function for each (p, phi) grid point
        for l = 1:gridSize
            for m = 1:gridSize
                norm_val = normpdf(0:N, N*p(l), sqrt(N*p(l)*(1-p(l))));
                phase = exp(-1i * (N:-1:0) * phi(m));
                HussimiMatrix(l, m) = abs(sqrt(norm_val) .* phase * v)^2;
            end
        end
        
        % Normalize the Husimi matrix for this state and add it to the sum
        Hussimi = Hussimi + HussimiMatrix / max(HussimiMatrix(:));
    end

    % Plot the accumulated and normalized Husimi function
    figure;
    hold on;
    imagesc([-1, 1], [-1, 1], real(Hussimi));
    set(gca, 'YDir', 'normal');
    pSpacePortrait(param(1), param(2), param(3));
    axis([-1 1 -1 1]);
    xlabel('$\varphi / \pi$', 'Interpreter', 'latex');
    ylabel('$p$', 'Interpreter', 'latex');
    hold off;
end

%%
% Function: Coherence
%
% Description:
%   Computes the coherence measure eta (as described in the publication) for a set 
% of Floquet states.
%
% Inputs:
%   - u_n: Matrix containing the Floquet states (each column is a state).
%   - N: Particle number.
%
function [eta, Indx] = Coherence(u_n, N)
    % Pre-calculate coefficients for one-particle reduced density matrix
    a11 = (0:N)'; 
    a22 = (N:-1:0)';
    a12 = (sqrt(N:-1:1) .* sqrt(1:N))';
    a21 = (sqrt(N:-1:1) .* sqrt(1:N))';
    
    for l = 1:N+1
        V = u_n(:, l);
        a1 = real(dot(V, a11 .* V));
        a2 = real(dot(V(1:end-1), a12 .* V(2:end)));
        a3 = real(dot(V(2:end), a21 .* V(1:end-1)));
        a4 = real(dot(V, a22 .* V));
        mu(l) = 2/N^2 * trace([a1, a2; a3, a4]^2) - 1;
    end
    [eta, Indx] = sort(mu);
end