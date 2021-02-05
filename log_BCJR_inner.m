% (inner decoder)log-BCJR algorithm
% outputs extrinsic information
% This log BCJR code is inspired by the BCJR (probability domain) scilab code written by K
% Vasudevan which can be found at http://home.iitk.ac.in/~vasu
function [LLR]= log_BCJR_inner(LLR,log_gamma,num_bit)
[Prev_State,Prev_Ip,Outputs_prev,Next_State,Outputs_next,Next_Ip]= Get_Trellis_Pred();
num_states =8; % number of states
C = exp(-1*LLR/2)./(1+exp(-1*LLR)); % Ck
repmat_C = repmat(C,num_states,1);
%******************************************************************************
% Initialize log-alpha and log-beta (assuming receiver does not know
% starting and ending states
%******************************************************************************
 log_alpha=zeros(num_states,num_bit);
 log_beta=zeros(num_states,num_bit+1);
 log_alpha(:,1)= 0;%  initialization
 log_beta(:,num_bit+1)= 0; % initialization
%******************************************************************************
%   Compute log-alpha and log-beta
%******************************************************************************
 for time=1:num_bit-1
     % forward recursion
     temp1 = log_alpha(Prev_State(:,1),time)+log(repmat_C(time))+(1-2*(Prev_Ip(:,1)-1))*LLR(time)/2+log_gamma(Outputs_prev(:,1),time);
     temp2 = log_alpha(Prev_State(:,2),time)+log(repmat_C(time))+(1-2*(Prev_Ip(:,2)-1))*LLR(time)/2+log_gamma(Outputs_prev(:,2),time);
     log_alpha(:,time+1)= max(temp1,temp2)+log(1+exp(-abs(temp1-temp2))) ; % Jacobian logarithm
     % backward recursion
     temp3 = log_beta(Next_State(:,1),num_bit+2-time)+log(repmat_C(num_bit+1-time))+(1-2*(Next_Ip(:,1)-1))*LLR(num_bit+1-time)/2+log_gamma(Outputs_next(:,1),num_bit+1-time);
     temp4 = log_beta(Next_State(:,2),num_bit+2-time)+log(repmat_C(num_bit+1-time))+(1-2*(Next_Ip(:,2)-1))*LLR(num_bit+1-time)/2+log_gamma(Outputs_next(:,2),num_bit+1-time);
     log_beta(:,num_bit+1-time)= max(temp3,temp4)+log(1+exp(-abs(temp3-temp4))) ; % Jacobian logarithm
 end

%**************************************************************************
% Compute extrinsic information
%**************************************************************************
 temp5 = log_alpha + log_gamma(Outputs_next(:,1),:)+ log_beta(Next_State(:,1),2:num_bit+1) ;
 temp5_1 = max(temp5(1,:),temp5(2,:))+log(1+exp(-abs(temp5(1,:)-temp5(2,:)))); % Jacobian logarithm
 temp5_2 = max(temp5(3,:),temp5(4,:))+log(1+exp(-abs(temp5(3,:)-temp5(4,:)))); % Jacobian logarithm
 temp5_3 = max(temp5(5,:),temp5(6,:))+log(1+exp(-abs(temp5(5,:)-temp5(6,:)))); % Jacobian logarithm
 temp5_4 = max(temp5(7,:),temp5(8,:))+log(1+exp(-abs(temp5(7,:)-temp5(8,:)))); % Jacobian logarithm
 temp5_11 = max(temp5_1,temp5_2)+log(1+exp(-abs(temp5_1-temp5_2))); % Jacobian logarithm
 temp5_22 = max(temp5_3,temp5_4)+log(1+exp(-abs(temp5_3-temp5_4))); % Jacobian logarithm
 LLR_temp1 = max(temp5_11,temp5_22)+log(1+exp(-abs(temp5_11-temp5_22))); % Jacobian logarithm
 
 temp6 = log_alpha + log_gamma(Outputs_next(:,2),:)+ log_beta(Next_State(:,2),2:num_bit+1) ;
 temp6_1 = max(temp6(1,:),temp6(2,:))+log(1+exp(-abs(temp6(1,:)-temp6(2,:)))); % Jacobian logarithm
 temp6_2 = max(temp6(3,:),temp6(4,:))+log(1+exp(-abs(temp6(3,:)-temp6(4,:)))); % Jacobian logarithm
 temp6_3 = max(temp6(5,:),temp6(6,:))+log(1+exp(-abs(temp6(5,:)-temp6(6,:)))); % Jacobian logarithm
 temp6_4 = max(temp6(7,:),temp6(8,:))+log(1+exp(-abs(temp6(7,:)-temp6(8,:)))); % Jacobian logarithm
 temp6_11 = max(temp6_1,temp6_2)+log(1+exp(-abs(temp6_1-temp6_2))); % Jacobian logarithm
 temp6_22 = max(temp6_3,temp6_4)+log(1+exp(-abs(temp6_3-temp6_4))); % Jacobian logarithm
 LLR_temp2 = max(temp6_11,temp6_22)+log(1+exp(-abs(temp6_11-temp6_22))); % Jacobian logarithm
 
 LLR = LLR_temp1 - LLR_temp2;

% normalizing to avoid numerical instabilities
LLR(LLR>50) = 50;
LLR(LLR<-50) = -50;
end