% Trellis for the Convolutional encoder with generator matrix given by
% [1 (1+D^2)/(1+D+D^2)]
function [Prev_State,Prev_Ip,Outputs_prev,Next_State,Outputs_next,Gb_Inx_par,Next_State_par]= Get_Trellis_Encoder()
Outputs_prev = [1,4; 2,3; 1,4; 2,3];  % Gamma indices for alpha recursion
Prev_State = [1,2; 4,3; 2,1; 3,4]; % previous state
Prev_Ip =  [1,2;1,2;1,2;1,2] ; % Previous input.
Next_State = [1,3; 3,1; 4,2; 2,4]; % Next state
Outputs_next = [1,4; 1,4; 2,3; 2,3]; % Gamma indices for beta recursion
Gb_Inx_par = [1,4;1,4;3,2;3,2]; % for parity bit
Next_State_par = [1,3;3,1;2,4;4,2]; % for parity bit
end
