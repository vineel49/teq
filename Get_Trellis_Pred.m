% Super Trellis for prediction filter order 1
% Generator matrix is G(D) = [1 (1+D^2)/(1+D+D^2)]
function [P_State,P_Ip,Ga_Inx,N_State,Gb_Inx,N_Ip]= Get_Trellis_Pred()
P_State = [1,2;3,4;7,8;5,6;3,4;1,2;5,6;7,8]; % previous state
P_Ip =  [1,1;2,2;1,1;2,2;1,1;2,2;1,1;2,2] ; % Previous input.
N_State = [1,6;1,6;5,2;5,2;7,4;7,4;3,8;3,8]; % Next state
Gb_Inx = [1,2;3,4;5,6;7,8;9,10;11,12;13,14;15,16]; % Gamma indices for beta recursion
Ga_Inx = [1,3;6,8;13,15;10,12;5,7;2,4;9,11;14,16];  % Gamma indices for alpha recursion
N_Ip = [1,2;1,2;1,2;1,2;1,2;1,2;1,2;1,2]; % Next input
end
