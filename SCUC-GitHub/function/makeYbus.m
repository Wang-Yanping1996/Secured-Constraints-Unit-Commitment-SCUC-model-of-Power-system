function [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch)
%MAKEYBUS   Builds the bus admittance matrix and branch admittance matrices.
%   [Ybus, Yf, Yt] = makeYbus(baseMVA, bus, branch) returns the full
%   bus admittance matrix (i.e. for all buses) and the matrices Yf and Yt
%   which, when multiplied by a complex voltage vector, yield the vector
%   currents injected into each line from the "from" and "to" buses
%   respectively of each line. Does appropriate conversions to p.u.

%   MATPOWER Version 2.0
%   by Ray Zimmerman, PSERC Cornell    12/19/97
%   Copyright (c) 1996, 1997 by Power System Engineering Research Center (PSERC)
%   See http://www.pserc.cornell.edu/ for more info.

%% constants
j = sqrt(-1);
nb = size(bus, 1);			%% number of buses
nl = size(branch, 1);		%% number of lines

%% ----- define named indices into bus, branch matrices
PQ=1; PV=2; REF=3; NONE=4; BUS_I=1; BUS_TYPE=2; PD=3; QD=4; GS=5; BS=6; BUS_AREA=7;
VM=8; VA=9; BASE_KV=10;	ZONE=11; VMAX=12; VMIN=13; LAM_P=14; LAM_Q=15; MU_VMAX=16; MU_VMIN=17; 

F_BUS=1; T_BUS=2; BR_R=3; BR_X=4; BR_B=5; RATE_A=6; RATE_B=7; RATE_C=8;% standard notation (in input)
TAP=9;  SHIFT=10; BR_STATUS=11; BR_angmin=12; BR_angmax=13;% standard notation (in input)
PF=12; QF=13; PT=14; QT=15; MU_SF=16; MU_ST=17; % MU_SF: idx of MU on MVA lim at "f" bus (in opf)

%% ----- check that bus numbers are equal to indices to bus (one set of bus numbers)
if any(bus(:, BUS_I) ~= [1:nb]')
	error('buses must appear in order by bus number')
end

%% for each branch, compute the elements of the branch admittance matrix where
%%
%%		| If |   | Yff  Yft |   | Vf |;
%%		|    | = |          | * |    |;
%%		| It |   | Ytf  Ytt |   | Vt |;
%%
stat = branch(:, BR_STATUS);					%% ones at in-service branches
Ys = stat ./ (branch(:, BR_R) + j * branch(:, BR_X));	%% series admittance
Bc = stat .* branch(:, BR_B);							%% line charging susceptance
tap = ones(nl, 1);								%% default tap ratio = 1
i = find(branch(:, TAP));						%% indices of non-zero tap ratios
tap(i) = branch(i, TAP);						%% assign non-zero tap ratios
tap = tap .* exp(j*pi/180 * branch(:, SHIFT));	%% add phase shifters
Ytt = Ys + j*Bc/2;
Yff = Ytt ./ (tap .* conj(tap));
Yft = - Ys ./ conj(tap);
Ytf = - Ys ./ tap;

%% compute shunt admittance
%% if Ps is the real power consumed by the shunt at V = 1.0 p.u.
%% and Qs is the reactive power injected by the shunt at V = 1.0 p.u.
%% then Ps - j Qs = V * conj(Ys * V) = conj(Ys) = Gs - j Bs,
%% i.e. Ys = Ps + j Qs, so ...
Ys = (bus(:, GS) + j * bus(:, BS)) / baseMVA;	%% vector of shunt admittances

%% build Ybus
f = branch(:, F_BUS);							%% list of "from" buses
t = branch(:, T_BUS);							%% list of "to" buses
Cf = sparse(f, 1:nl, ones(nl, 1), nb, nl);		%% connection matrix for line & from buses
Ct = sparse(t, 1:nl, ones(nl, 1), nb, nl);		%% connection matrix for line & to buses
Ybus = spdiags(Ys, 0, nb, nb) + ...				%% shunt admittance
	Cf * spdiags(Yff, 0, nl, nl) * Cf' + ...	%% Yff term of branch admittance
	Cf * spdiags(Yft, 0, nl, nl) * Ct' + ...	%% Yft term of branch admittance
	Ct * spdiags(Ytf, 0, nl, nl) * Cf' + ...	%% Ytf term of branch admittance
	Ct * spdiags(Ytt, 0, nl, nl) * Ct';			%% Ytt term of branch admittance

%% Build Yf and Yt such that Yf * V is the vector of complex branch currents injected
%% at each branch's "from" bus, and Yt is the same for the "to" bus end
if nargout > 1
	i = [[1:nl]'; [1:nl]'];		%% double set of row indices	
	Yf = sparse(i, [f; t], [Yff; Yft]);
	Yt = sparse(i, [f; t], [Ytf; Ytt]);
end

return;

