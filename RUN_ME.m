% Turbo  equalization  of  serially  concatenated  turbo  codes using a predictive DFE-based receiver

close all
clear all
clc
%---------------- SIMULATION PARAMETERS ------------------------------------
SNR_dB = 8; % SNR per bit in dB (in logarithmic scale)
sim_runs = 1*(10^1); % simulation runs
training_len = 1000; % length of the training sequence
equalizer_len = 30; % length of the LMMSE equalizer
frame_size = 1000; % frame size
num_bit = 0.5*frame_size; % number of data bits (overall rate is 1/2)
fade_chan = [0.9 0.1 0.1]+1i*[0.9 0.1 0.1]; % ISI CHANNEL
fade_chan = fade_chan/norm(fade_chan); % normalizing to unit energy
SNR = 10^(0.1*SNR_dB); % SNR per bit in linear scale
noise_var_1D = 2*2/(2*SNR); % 1D noise variance
%--------------------------------------------------------------------------

%    Generator polynomial of the inner encoder
gen_poly_inner = ldiv2([1 0 1],[1 1 1],2*num_bit); % using long division method

%    Generator polynomial of the outer encoder
gen_poly_outer = ldiv2([1 0 1],[1 1 1],num_bit); % using long division method

%  Interleaver and deinterleaver mapping of the SCCC 
intr_map = randperm(2*num_bit);
deintr_map = deintrlv((1:2*num_bit),intr_map);
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                       TRAINING PHASE
% Source
training_a = randi([0 1],1,2*training_len); % data

% QPSK mapping 
training_sig = 1-2*training_a(1:2:end) + 1i*(1-2*training_a(2:2:end));

%                            CHANNEL   
% AWGN
white_noise = sqrt(noise_var_1D)*randn(1,training_len+length(fade_chan)-1)+1i*sqrt(noise_var_1D)*randn(1,training_len+length(fade_chan)-1); 
chan_op = conv(training_sig,fade_chan) + white_noise; % equalizer_op stands for channel output

%                           RECEIVER (TRAINING)
% Estimating the autocorrelation of equalizer_op at zero lag. 
Rvv0 = (chan_op*chan_op')/(training_len+length(fade_chan)-1);

%    LMS update of taps
equalizer = zeros(1,equalizer_len);
max_step_size = 2/(equalizer_len*Rvv0);% maximum step size
step_size = 0.5*max_step_size; % step size

for i1=1:training_len-equalizer_len+1
    equalizer_ip = fliplr(chan_op(i1:i1+equalizer_len-1));%equalizer input
    error = training_sig(i1+equalizer_len-1)- equalizer*equalizer_ip.';% instantaneous error
    equalizer = equalizer + step_size*error*conj(equalizer_ip);
end

% estimating the autocorrelation of the noise at the equalizer output
temp = conv(equalizer,chan_op);
equalizer_op = temp(1:training_len); % equalizer output

corr_vec = zeros(3,1); % autocorrelation vector initialization
noise_eq_op = equalizer_op-training_sig; % noise at equalizer output
corr_vec(1) = real(noise_eq_op*noise_eq_op'/(2*training_len)); % at zero lag
corr_vec(2) = noise_eq_op(2:end)*noise_eq_op(1:end-1)'/(2*(training_len-1));
corr_vec(3) = noise_eq_op(3:end)*noise_eq_op(1:end-2)'/(2*(training_len-2));

% Generate linear prediction filter coefficients using Levinson-Durbin
% algorithm
[pred_coef,pred_var,ref_coef]=Gen_Coef(corr_vec,1); % first order

%--------------------------------------------------------------------------
C_Ber = 0; % channel erros
tic()
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%                            DATA PHASE
for frame_cnt = 1:sim_runs
%                           TRANSMITTER
%Source
a = randi([0 1],1,num_bit); % data

% SCCC encoder
% Outer encoder
b = zeros(1,2*num_bit); % outer encoder output initialization
b(1:2:end) = a; % systematic bit
temp1 = mod(conv(gen_poly_outer,a),2); % linear convolution with the generator polynomial
b(2:2:end) = temp1(1:num_bit); % parity bit

% interleaver
c = b(intr_map);

% Inner encoder
d = zeros(1,2*frame_size); % inner encoder output initialization
d(1:2:end) = c; % systematic bit
temp2 = mod(conv(gen_poly_inner,c),2); % linear convolution with the generator polynomial
d(2:2:end) = temp2(1:frame_size); % parity bit

% QPSK mapping (according to the set partitioning principles)
mod_sig = 1-2*d(1:2:end) + 1i*(1-2*d(2:2:end));

%--------------------------------------------------------------------------
%                            CHANNEL   
% AWGN
white_noise = sqrt(noise_var_1D)*randn(1,frame_size+length(fade_chan)-1)+1i*sqrt(noise_var_1D)*randn(1,frame_size+length(fade_chan)-1); 
chan_op = conv(mod_sig,fade_chan) + white_noise; % equalizer_op stands for channel output
%--------------------------------------------------------------------------
%                          RECEIVER 
% LMMSE EQUALIZER
temp = conv(chan_op,equalizer);
equalizer_op = temp(1:frame_size);

% Branch metrices for the inner BCJR
QPSK_SYM = [1+1i 1-1i -1+1i -1-1i];
Dist = zeros(16,frame_size);
Dist(1,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(1))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(1)))).^2;
Dist(2,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(4))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(1)))).^2;

Dist(3,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(1))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(4)))).^2;
Dist(4,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(4))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(4)))).^2;

Dist(5,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(1))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(2)))).^2;
Dist(6,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(4))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(2)))).^2;

Dist(7,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(1))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(3)))).^2;
Dist(8,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(4))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(3)))).^2;

Dist(9,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(2))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(1)))).^2;
Dist(10,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(3))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(1)))).^2;

Dist(11,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(2))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(4)))).^2;
Dist(12,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(3))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(4)))).^2;

Dist(13,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(2))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(2)))).^2;
Dist(14,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(3))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(2)))).^2;

Dist(15,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(2))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(3)))).^2;
Dist(16,2:end) = abs((equalizer_op(2:end)-QPSK_SYM(3))+(pred_coef*(equalizer_op(1:end-1)-QPSK_SYM(3)))).^2;
log_gamma = -Dist/(2*pred_var); % log gamma
 
% a priori LLR for inner decoder for 1st iteration
LLR = zeros(1,frame_size);

% iterative logMAP decoding
LLR = log_BCJR_inner(LLR,log_gamma,frame_size); % outputs extrinsic information
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %1

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %2

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %3

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %4

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %5

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %6

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size); 
LLR = log_BCJR_outer(LLR(deintr_map),num_bit); %7

LLR = log_BCJR_inner(LLR(intr_map),log_gamma,frame_size);
LLR = log_BCJR_outer_END(LLR(deintr_map),num_bit); % 8: outputs aposteriori information

% hard decision 
dec_data = LLR<0;
% 
 % Calculating total bit errors
C_Ber = C_Ber + nnz(dec_data-a); 
end

BER = C_Ber/(sim_runs*num_bit)
toc()

