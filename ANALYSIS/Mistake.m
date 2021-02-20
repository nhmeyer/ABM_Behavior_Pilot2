function [hitcount,misscount,Correct_T_Env1,Wrong_T_Env1] = Mistake(NumQ1,IndexT_Env1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

a = 1; b = 1; hitcount = 0; misscount = 0;
for i = 1:length(IndexT_Env1)
    if NumQ1(IndexT_Env1(i),1) == 0 
        hitcount = hitcount + 1
        Hit_T = hitcount;
        Correct_T_Env1(a) = i;
        a = a+1;
    else
        misscount = misscount+1
        Miss_T = misscount
        Wrong_T_Env1(b) = i;
        
        b = b+1;
       
        end
end
if ~exist('Wrong_T_Env1','var')
    Wrong_T_Env1 =0;
end
if ~exist('Correct_T_Env1','var')
   Correct_T_Env1 =0;
end
end

