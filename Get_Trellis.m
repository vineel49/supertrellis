% differential detection - for convolutional coded DQPSK in SIMO-OFDM
% prediction order 3, 16 states.

function[Next_State,Prev_State,Prev_State_trans,Prev_Ip,Prev_Ip_trans,Outputs_next,Outputs_prev] = Get_Trellis()

Next_State = [1,11; 1,11; 2,12; 2,12; 9,3; 9,3; 10,4; 10,4; 13,7; 13, 7; 14,8; 14,8; 5, 15; 5, 15; 6, 16; 6, 16];

Prev_State = [1,2; 3,4; 5,6; 7,8; 13,14; 15,16; 9,10; 11,12; 5,6; 7,8; 1,2; 3,4; 9,10; 11, 12; 13,14; 15,16];

Prev_State_trans = [1,3,5,7,13,15,9,11,5,7,1,3,9,11,13,15;2,4,6,8,14,16,10,12,6,8,2,4,10,12,14,16];

Prev_Ip = [1,1; 1,1; 2,2; 2,2; 1,1; 1,1; 2,2; 2,2; 1,1; 1,1; 2,2; 2,2; 1,1; 1,1; 2,2; 2,2];

Prev_Ip_trans = [1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2;1,1,2,2,1,1,2,2,1,1,2,2,1,1,2,2];

Outputs_next = [1,2;3,4;5,6;7,8;9,10;11,12;13,14;15,16;17,18;19,20;21,22;23,24;25,26;27,28;29,30;31,32];

Outputs_prev = [1,3;5,7;10,12;14,16;25,27;29,31;18,20;22,24;9,11;13,15;2,4;6,8;17,19;21,23;26,28;30,32]; 
end