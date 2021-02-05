% Super Trellis for prediction filter order 1
% Generator matrix is G(D) = [1 (1+D^2)/(1+D+D^2)]
function [Prev_State,Prev_Ip,Outputs_prev,Next_State,Outputs_next,Next_Ip]= Get_Trellis_Pred()
Prev_State = [1,2;3,4;7,8;5,6;3,4;1,2;5,6;7,8]; % previous state
Prev_Ip =  [1,1;2,2;1,1;2,2;1,1;2,2;1,1;2,2] ; % Previous input.
Next_State = [1,6;1,6;5,2;5,2;7,4;7,4;3,8;3,8]; % Next state
Outputs_next = [1,2;3,4;5,6;7,8;9,10;11,12;13,14;15,16]; % Gamma indices for beta recursion
Outputs_prev = [1,3;6,8;13,15;10,12;5,7;2,4;9,11;14,16];  % Gamma indices for alpha recursion
Next_Ip = [1,2;1,2;1,2;1,2;1,2;1,2;1,2;1,2]; % Next input
end
