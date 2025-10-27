function [V, delta, Psl, Qgv, N, time] = oliver_nrpf(Y, is, ipq, ipv, Pg, Qg, Pd, Qd, V0, Sbase, toler, maxiter)

% setup Y, etc

tic;

d0 = 0;

N = length(Y);

jp = setdiff(1:N, is)';

jq = ipq;


d = zeros(N, 1);   % delta

Vm = zeros(N, 1);  % voltage magnitudes


% initialize values with flat-voltage profile


d(is) = d0;

Vm(ipv) = V0; 
Vm(jq) = Vm(jq) + 1;

V = Vm .* exp(j*d);


Pinj = (Pg - Pd) / Sbase;

Qinj = (Qg - Qd) / Sbase;


S = V .* conj(Y * V);

change = [Pinj(jp) - real(S(jp)); Qinj(jq) - imag(S(jq))];


k = 1;

while k <= maxiter && max(abs(change)) >= toler
    
    % calculate jacobian and change

    J1 = real(-j*diag(V)*conj(Y)*diag(conj(V)) + j*diag(conj(Y)*conj(V))*diag(V));
    J2 = real(diag(V)*conj(Y)*diag(exp(-j*d)) + diag(conj(Y)*conj(V))*diag(exp(j*d)));
    J3 = imag(-j*diag(V)*conj(Y)*diag(conj(V)) + j*diag(conj(Y)*conj(V))*diag(V));
    J4 = imag(diag(V)*conj(Y)*diag(exp(-j*d)) + diag(conj(Y)*conj(V))*diag(exp(j*d)));

    J = [J1(jp,jp), J2(jp,jq); J3(jq,jp), J4(jq,jq)];


    xchange = J \ change;


    % update values

    d(jp) = d(jp) + xchange(1:length(jp));
    Vm(jq) = Vm(jq) + xchange(length(jp)+1:end);
    V = Vm .* exp(j*d);

    S = V .* conj(Y * V);

    change = [Pinj(jp) - real(S(jp)); Qinj(jq) - imag(S(jq))];

    
    k = k + 1;

end

time = toc;

if k > maxiter + 1
    error('Power flow calculation did not converge within the maximum number of iterations.');

else
    
    N = k - 1;
    
    V = Vm;

    delta = d;

    Psl = Pinj(1) + real(S(1));

    Qgv = Qinj + imag(S);


end    


end