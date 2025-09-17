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