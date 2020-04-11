% Viterbi Algorithm - Soft Input, Hard Output 
% Written by Vineel Kumar Veludandi
function[dec_ip]=Viterbi_algorithm(num_sym,decoding_delay,branch_metric)
[~,P_State,P_State_trans,P_Ip,~,~,Ga_Inx] = Get_Trellis();
num_states = 16;  % number of states in the convolutional encoder
% General Initialization
ip=0; % input initialization
dec_ip = zeros(1,num_sym-decoding_delay); % decoded symbols (-decoding delay is because,ignoring last transient part)
survivor_node = zeros(num_states,num_sym); % survivor nodes
survivor_ip = zeros(num_states,num_sym); % survivor inputs
path_metric = zeros(num_states,num_sym+1); % path metrics
index_temp = [0;1*2;2*2;3*2;4*2;5*2;6*2;7*2;8*2;9*2;10*2;11*2;12*2;13*2;14*2;15*2]; %for linear indexing.
for sym_cnt=  1:num_sym  
   [path_metric(:,sym_cnt+1),index] = min([path_metric(P_State(:,1),sym_cnt)+ branch_metric(Ga_Inx(:,1),sym_cnt) ...
       path_metric(P_State(:,2),sym_cnt)+ branch_metric(Ga_Inx(:,2),sym_cnt)],[],2);
   survivor_node(:,sym_cnt) = P_State_trans(index+index_temp);
   survivor_ip(:,sym_cnt) = P_Ip(index+index_temp); 
   %      Back tracing  
if (sym_cnt>decoding_delay) 
[~,trace_bf] = min(path_metric(:,sym_cnt+1));
for bk_cnt= 1 : decoding_delay
ip = survivor_ip(trace_bf,sym_cnt+1-bk_cnt);
trace_bf = survivor_node(trace_bf,sym_cnt+1-bk_cnt);    
end
dec_ip(sym_cnt-decoding_delay)=ip;
end % for if statement
end % for forward recursiion
dec_ip = dec_ip-1; % decoded bits ( now in 0 and 1)
end % for function
