function [F1_r,F2_r,F3_r,F4_r,F5_r,F6_r,F7_r,F8_r] = MistakeDetailed(myexcelfile,MyParticipants,i,myEnv)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[Num,Text,All] = xlsread(myexcelfile)
Text_1 = contains(Text(:,1),'Q1]');
 NumQ1 = Num(find(Text_1),:)
 TextQ1 = Text(Text_1,:);
IndexEnv1 = contains(TextQ1(:,2),sprintf('%s',myEnv));nnz(IndexEnv1)
NumQ1Env1 = NumQ1(IndexEnv1,1);
F1 = contains(TextQ1(IndexEnv1,2),'F1');
F2 = contains(TextQ1(IndexEnv1,2),'F2');
F3 = contains(TextQ1(IndexEnv1,2),'F3');
F4 = contains(TextQ1(IndexEnv1,2),'F4');
F5 = contains(TextQ1(IndexEnv1,2),'F5')
F6 = contains(TextQ1(IndexEnv1,2),'F6');
F7 = contains(TextQ1(IndexEnv1,2),'F7');
F8 = contains(TextQ1(IndexEnv1,2),'F8');
NumF1 = NumQ1Env1((F1));
NumF2 = NumQ1Env1((F2));
NumF3 = NumQ1Env1((F3));
NumF4 = NumQ1Env1((F4));
NumF5 = NumQ1Env1((F5));
NumF6 = NumQ1Env1((F6));
NumF7 = NumQ1Env1((F7));
NumF8 = NumQ1Env1((F8));
a = 0;
for i = 1:length(NumF1)
    if NumF1(i) == 0
    a =a+1;
    end
     MyF1 = a;
end
a = 0;
for i = 1:length(NumF7)
    if NumF7(i) == 0
         a =a+1;
    end
     MyF7 = a;
end

a = 0;
for i = 1:length(NumF6)
    if NumF6(i) == 0
          a =a+1;
    end
     MyF6 = a;
end

a = 0;
for i = 1:length(NumF5)
    if NumF5(i) == 0
          a =a+1;
    end
     MyF5 = a;
end

a = 0;
for i = 1:length(NumF4)
    if NumF4(i) == 0
          a =a+1;
    end
     MyF4 = a;
end

a = 0;
for i = 1:length(NumF3)
    if NumF3(i) == 0
           a =a+1;
    end
     MyF3 = a;
end

a = 0;
for i = 1:length(NumF2)
    if NumF2(i) == 0
            a =a+1;
    end
     MyF2 = a;
end

a = 0;
for i = 1:length(NumF8)
    if NumF8(i) == 0
          a =a+1;
    end
     MyF8 = a;
end

a = 0;
if isempty(NumF2)
    NumF2 = 1;
    MyF2 = 0;
end
F1_r = MyF1/length(NumF1);F2_r = MyF2/length(NumF2);F3_r = MyF3/length(NumF3);F4_r = MyF4/length(NumF4);F5_r = MyF5/length(NumF5);F6_r = MyF6/length(NumF6);F7_r = MyF7/length(NumF7);F8_r = MyF8/length(NumF8);
end

