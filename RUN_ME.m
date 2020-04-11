% Linear Prediction based Data Detection of Convolutional Coded 
% DQPSK in SIMO-OFDM
% Prediction order 3, 16 states. (convolutional coded DQPSK modulation)
% for Diversity orders 4 , frame size 1024. Overall rate is 1.
% Written by Vineel Kumar Veludandi

close all
clear all
clc
SNR_dB = 20; % SNR per bit in dB (in logarithmic scale)
frame_size = 1024; % 1024 pt. FFT and IFFT length
chan_len = 10; % number of channel taps
num_frames = 1*(10^1); % number of frames to be simulated. - frame count
cp_num = chan_len-1; % length of the cyclic prefix
fade_var_1D =0.5; % 1D fade variance
fade_std_dev =sqrt(fade_var_1D);
num_bit = frame_size*1; % number of data bits, overall rate 1
decoding_delay = 20; % decoding delay of the Viterbi algorithm

% Computing impulse response of the RSC encoder using long division method
gen_poly = ldiv2([1 0 1],[1 1 1],num_bit); 

% SNR parameters - overall rate 1 - receiver diversity order is 4
SNR = 10^(0.1*SNR_dB); % converting logarithmic scale into linear scale
NOISE_VAR_1D = 4*chan_len/(frame_size*SNR); % 1D noise variance
NOISE_STD_DEV = sqrt(NOISE_VAR_1D);

% Generate the channel DFT autocorrelation from analytical expression (see shimla paper)
auto_cor = zeros(6,1);
for auto_cnt = 0:chan_len-1
    auto_cor(1) = auto_cor(1) + (fade_var_1D *exp(-1i*2*pi*0*auto_cnt/frame_size));
    auto_cor(2) = auto_cor(2) + (fade_var_1D *exp(-1i*2*pi*1*auto_cnt/frame_size));
    auto_cor(3) = auto_cor(3) + (fade_var_1D *exp(-1i*2*pi*2*auto_cnt/frame_size));
    auto_cor(4) = auto_cor(4) + (fade_var_1D *exp(-1i*2*pi*3*auto_cnt/frame_size));
    auto_cor(5) = auto_cor(5) + (fade_var_1D *exp(-1i*2*pi*4*auto_cnt/frame_size));
    auto_cor(6) = auto_cor(6) + (fade_var_1D *exp(-1i*2*pi*5*auto_cnt/frame_size));
end
    
% taking noise variance into account
auto_cor(1) = auto_cor(1) + (NOISE_VAR_1D *frame_size/2);
%-----------------------------------------------------------------------------

% generate prediction coefficients based on auto correlation alone
[pred_coef,pred_error] = Gen_Coef(auto_cor,3);
[pred_coef_2tap,pred_error_2tap] = Gen_Coef(auto_cor,2);
[pred_coef_1tap,pred_error_1tap] = Gen_Coef(auto_cor,1);

C_Ber=0; % initialization for bit errors - total
tic()
for frame_cnt = 1:num_frames
   
% source bit generation
a= randi([0 1],1,num_bit);

% CONVOLUTIONAL ENCODER
b = zeros(1,2*num_bit); % ENCODER OUTPUT INITIALIZATION
b(1:2:end) = a; % SYSTEMATIC BIT
temp = mod(conv(gen_poly,a),2); 
b(2:2:end) = temp(1:num_bit); % PARITY BIT 

% dqpsk mapping
dqpsk_seq = zeros(1,num_bit);
ref_sym = 1+1i; % reference symbol
dqpsk_rules = [exp(1i*0) exp(1i*pi/2) exp(1i*3*pi/2) exp(1i*pi)]; % dqpsk encoding rules
dqpsk_seq(1) = ref_sym * dqpsk_rules(2*b(1)+b(2)+1);
for i = 2:(num_bit) 
   dqpsk_seq(i) = dqpsk_seq(i-1)* dqpsk_rules( 2*b(2*i-1)+b(2*i)+1);
end   

%--------------------------------------------------------------------------
% OFDM Transmitter
F_trans_sig_no_CP = dqpsk_seq; % F stands for frequency domain

% taking IFFT 
T_trans_sig_no_CP = ifft(F_trans_sig_no_CP); % T for time domain - IFFT

% adding CYCLIC PREFIX
T_trans_sig_CP = [T_trans_sig_no_CP(end - cp_num+1:end) T_trans_sig_no_CP]; % adding cyclic prefix. 
%-------------------------------------------------------------------------- 
% Get ISI channel
fade_mean = 0;
fade_chan1 = normrnd(0,fade_std_dev,1,chan_len)+1i*normrnd(0,fade_std_dev,1,chan_len);     
fade_chan2 = normrnd(0,fade_std_dev,1,chan_len)+1i*normrnd(0,fade_std_dev,1,chan_len);     
fade_chan3 = normrnd(0,fade_std_dev,1,chan_len)+1i*normrnd(0,fade_std_dev,1,chan_len);     
fade_chan4 = normrnd(0,fade_std_dev,1,chan_len)+1i*normrnd(0,fade_std_dev,1,chan_len);   
    
%--------------------------------------------------------------------------
Chan_Op_temp1 = conv(T_trans_sig_CP,fade_chan1); % convolution with the channel    
Chan_Op_temp2 = conv(T_trans_sig_CP,fade_chan2); % convolution with the channel    
Chan_Op_temp3 = conv(T_trans_sig_CP,fade_chan3); % convolution with the channel    
Chan_Op_temp4 = conv(T_trans_sig_CP,fade_chan4); % convolution with the channel 
    
% AWGN noise generation
noise1 = normrnd(0,NOISE_STD_DEV,1,frame_size+cp_num+chan_len-1)+1i*normrnd(0,NOISE_STD_DEV,1,frame_size+cp_num+chan_len-1);    
noise2 = normrnd(0,NOISE_STD_DEV,1,frame_size+cp_num+chan_len-1)+1i*normrnd(0,NOISE_STD_DEV,1,frame_size+cp_num+chan_len-1);    
noise3 = normrnd(0,NOISE_STD_DEV,1,frame_size+cp_num+chan_len-1)+1i*normrnd(0,NOISE_STD_DEV,1,frame_size+cp_num+chan_len-1);    
noise4 = normrnd(0,NOISE_STD_DEV,1,frame_size+cp_num+chan_len-1)+1i*normrnd(0,NOISE_STD_DEV,1,frame_size+cp_num+chan_len-1);  

T_rec_sig_CP1 = Chan_Op_temp1 + noise1;    
T_rec_sig_CP2 = Chan_Op_temp2 + noise2;    
T_rec_sig_CP3 = Chan_Op_temp3 + noise3;    
T_rec_sig_CP4 = Chan_Op_temp4 + noise4;  
%--------------------------------------------------------------------------
% Receiver Operations
% CP & transient samples removal 
T_rec_sig_CP1(1:cp_num) = [];
T_rec_sig_no_CP1 = T_rec_sig_CP1(1:frame_size); % cyclic prefix removal    

T_rec_sig_CP2(1:cp_num) = [];
T_rec_sig_no_CP2 = T_rec_sig_CP2(1:frame_size); % cyclic prefix removal    

T_rec_sig_CP3(1:cp_num) = [];
T_rec_sig_no_CP3 = T_rec_sig_CP3(1:frame_size); % cyclic prefix removal    

T_rec_sig_CP4(1:cp_num) = [];
T_rec_sig_no_CP4 = T_rec_sig_CP4(1:frame_size); % cyclic prefix removal   

    
%--------------------------------------------------------------------------
% performing the FFT
F_rec_sig_no_CP1 = fft(T_rec_sig_no_CP1); % performing the FFT operation
F_rec_sig_no_CP2 = fft(T_rec_sig_no_CP2); % performing the FFT operation
F_rec_sig_no_CP3 = fft(T_rec_sig_no_CP3); % performing the FFT operation
F_rec_sig_no_CP4 = fft(T_rec_sig_no_CP4); % performing the FFT operation

% Gamma generation
branch_metric1 = Gen_Gamma_Pred(frame_size,F_rec_sig_no_CP1,pred_coef_1tap,pred_coef_2tap,pred_coef);
branch_metric2 = Gen_Gamma_Pred(frame_size,F_rec_sig_no_CP2,pred_coef_1tap,pred_coef_2tap,pred_coef);
branch_metric3 = Gen_Gamma_Pred(frame_size,F_rec_sig_no_CP3,pred_coef_1tap,pred_coef_2tap,pred_coef);
branch_metric4 = Gen_Gamma_Pred(frame_size,F_rec_sig_no_CP4,pred_coef_1tap,pred_coef_2tap,pred_coef);
branch_metric = branch_metric1 + branch_metric2 +branch_metric3 + branch_metric4;

dec_a=Viterbi_algorithm(frame_size,decoding_delay,branch_metric);

% Considering only the steady-state part in the decision rule, hence considering a(3:$)- i.e., ignoring first dqpsk symbol
C_Ber = C_Ber + nnz((dec_a(3:end))-a(3:end-decoding_delay)); % ignoring transient state, thats why a(3:$)
end
toc()
BER = (C_Ber)/((num_bit-2-decoding_delay)*num_frames) % ignoring transient state
    




