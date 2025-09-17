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
    set(gca,'ydir','normal','fontsize',24);
    xlabel('$\varphi / \pi$', 'Interpreter', 'latex',FontSize=28);
    ylabel('$p$', 'Interpreter', 'latex',FontSize=28);
    hold off;
end