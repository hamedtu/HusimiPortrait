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