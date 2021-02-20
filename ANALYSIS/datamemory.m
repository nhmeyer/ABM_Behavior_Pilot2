MyParticipants = ls('DesktopVRComputeranalysis/P0*');
%MyParticipants = MyParticipants([2:9 11:end],:);
list_env1issync = [2 2 2 1 1 1 1 2 2 2 2 2 2 1 1 1 1 1 2 1 2 1];
for i = 2%:size(MyParticipants,1)
    
[Num,Text,All] = xlsread(sprintf('DesktopVRComputeranalysis/%s/%s',MyParticipants(i,:),ls(sprintf('DesktopVRComputeranalysis/%s/ROOMQ*.xlsx',MyParticipants(i,:)))))
[NumBSC,TextBSC,AllBSC] = xlsread(sprintf('DesktopVRComputeranalysis/%s/%s',MyParticipants(i,:),ls(sprintf('DesktopVRComputeranalysis/%s/BSCQ*.xlsx',MyParticipants(i,:)))))

Text_1 = contains(Text(:,1),'Q1]');nnz(Text_1)
TextQ1 = Text(Text_1,:);
Text_2 = contains(Text(:,1),'Q2]');nnz(Text_2)
TextQ2 = Text(Text_2,:);
 NumQ1 = Num(find(Text_1),:)
 NumQ2 = Num(find(Text_2),:)
IndexEnv1 = contains(TextQ1(:,2),'Env1]');nnz(IndexEnv1)
IndexEnv2 = contains(TextQ1(:,2),'Env2]');nnz(IndexEnv2)
NumQ1Env1 = NumQ1(IndexEnv1,:);
NumQ1Env2 = NumQ1(IndexEnv2,:);

IndexF_Env1 = contains(TextQ1(IndexEnv1,2),'F');
Index_F1_ENV1 = contains(TextQ1(IndexEnv1,2),'F1'); 
Index_F2_ENV1 = contains(TextQ1(IndexEnv1,2),'F2'); 

%IndexF_Env1([find(Index_F1_ENV1)' find(Index_F2_ENV1)']) = [] uncomment if
%you want to analyze this without the object F1 and F2 of the kitchen
%environement (statistiacly different in term of difficulty compared to the
%other objects)

Index_Catch_ENV1 = contains(TextQ1(IndexEnv1,2),'Catch'); 
Index_Catch_ENV2 = contains(TextQ1(IndexEnv2,2),'Catch');
Index_Catch_ENV2_2 = contains(TextQ1(IndexEnv2,2),'Pos'); % This trick because some  catch trials are not labeled with "Catch" but "Pos". In this line we find the index of the "Pos" labeled catchtrials 
Index_Catch_ENV2(find(Index_Catch_ENV2_2)) = 1; % In this line we report the index of the "Pos" catch trial on the other vector in order to have one vectore containing all the label of the Catch Trial.

TotFalse_Env1 = nnz(IndexF_Env1);
IndexT_Env1 = contains(TextQ1(IndexEnv1,2),'T');
TotTrue_Env1 = nnz(IndexT_Env1)
T_Env1 = find(IndexT_Env1);
F_Env1 = find(IndexF_Env1);
TotCatch_Env1 = nnz(Index_Catch_ENV1);
TotCatch_Env2 = nnz(Index_Catch_ENV2);
CATCH_ENV1 = find(Index_Catch_ENV1);
CATCH_ENV2 = find(Index_Catch_ENV2);



IndexF_Env2 = contains(TextQ1(IndexEnv2,2),'F');nnz(IndexF_Env2)
IndexT_Env2= contains(TextQ1(IndexEnv2,2),'T');nnz(IndexT_Env2)
TotTrue_Env2 = nnz(IndexT_Env2)
TotFalse_Env2 = nnz(IndexF_Env2);

T_Env2 = find(IndexT_Env2);
F_Env2 = find(IndexF_Env2);

% performance
[correctrej_Env1,falsealarm_Env1,Correct_T_Env1,Wrong_T_Env1] =Mistake(NumQ1Env1,T_Env1);
[correctrej_Env2,falsealarm_Env2,Correct_T_Env2,Wrong_T_Env2] =Mistake(NumQ1Env2,T_Env2);
[misscount_Env1,hitcount_Env1,Wrong_F_Env1,Correct_F_Env1] =Mistake(NumQ1Env1,F_Env1);
[misscount_Env2,hitcount_Env2,Wrong_F_Env2,Correct_F_Env2] =Mistake(NumQ1Env2,F_Env2);
[misscountCATCH_Env1,hitcountCATCH_Env1,~,~] =Mistake(NumQ1Env1,CATCH_ENV1);
[misscountCATCH_Env2,hitcountCATCH_Env2,~,~] =Mistake(NumQ1Env2,CATCH_ENV2);

[hitcount_Env1 hitcount_Env2; misscount_Env1 misscount_Env2;falsealarm_Env1 falsealarm_Env2; correctrej_Env1 correctrej_Env2]
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{'ParticipantsName'},{'HitEnv1'},{'MissEnv1'},{'FalseAlarmEnv1'},{'CorrectRejectionEnv1'},{'N_TrueEnv1'},{'N_FalseEnv1'},{'HitEnv2'},{'MissEnv2'},{'FalseAlarmEnv2'},{'CorrectRejectionEnv2'},{'N_TrueEnv2'},{'N_FalseEnv2'}],'Performance','A1:M1')
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',{MyParticipants(i,:)},'Performance',sprintf('A%s',num2str(i+1)))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[hitcount_Env1 misscount_Env1 falsealarm_Env1 correctrej_Env1 length(T_Env1) length(F_Env1) hitcount_Env2 misscount_Env2 falsealarm_Env2 correctrej_Env2 length(T_Env2) length(F_Env2)],'Performance',sprintf('B%s:M%s',num2str(i+1),num2str(i+1)))

%With CatchTrials
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{'ParticipantsName'},{'HitEnv1'},{'MissEnv1'},{'FalseAlarmEnv1'},{'CorrectRejectionEnv1'},{'N_TrueEnv1'},{'N_FalseEnv1'},{'HitEnv2'},{'MissEnv2'},{'FalseAlarmEnv2'},{'CorrectRejectionEnv2'},{'N_TrueEnv2'},{'N_FalseEnv2'}],'PerformanceCatch','A1:M1')
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',{MyParticipants(i,:)},'PerformanceCatch',sprintf('A%s',num2str(i+1)))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[hitcountCATCH_Env1 misscountCATCH_Env1 falsealarm_Env1 correctrej_Env1 length(T_Env1) TotCatch_Env1 hitcountCATCH_Env2 misscountCATCH_Env2 falsealarm_Env2 correctrej_Env2 length(T_Env2) TotCatch_Env2],'PerformanceCatch',sprintf('B%s:M%s',num2str(i+1),num2str(i+1)))


% Confidence
MeanConfidenceT_Env1 = nanmean(NumQ2(T_Env1,1));
MeanConfidenceT_Env2 = nanmean(NumQ2(T_Env2,1));
MeanConfidenceF_Env1 = nanmean(NumQ2(F_Env1,1));
MeanConfidenceF_Env2 = nanmean(NumQ2(F_Env2,1));
[MeanConfidenceT_Env1 MeanConfidenceT_Env2; MeanConfidenceF_Env1 MeanConfidenceF_Env2]
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{'ParticipantsName'},{'Confidence_T_Env1'},{'Confidence_F_Env1'},{'Confidence_T_Env2'},{'Confidence_F_Env2'}],'Confidence','A1:E1')
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',{MyParticipants(i,:)},'Confidence',sprintf('A%s',num2str(i+1)))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[MeanConfidenceT_Env1 MeanConfidenceF_Env1 MeanConfidenceT_Env2 MeanConfidenceF_Env2],'Confidence',sprintf('B%s:E%s',num2str(i+1),num2str(i+1)))
%  
%BSCQ
TextBSC_1 = contains(TextBSC(:,2),'ENV1');nnz(TextBSC_1)
TextBSC_Env1 = TextBSC(TextBSC_1,:);
TextBSC_2 = contains(TextBSC(:,2),'ENV2');nnz(TextBSC_2)
TextBSC_Env2 = TextBSC(TextBSC_2,:);
NumBSC_1 = NumBSC(find(TextBSC_1),:)
NumBSC_2 = NumBSC(find(TextBSC_2),:)
Ownership_1 = NumBSC_1(contains(TextBSC_Env1(:,2),'OWNERSHIP]'));nnz(Ownership_1)
Agency_1 = NumBSC_1(contains(TextBSC_Env1(:,2),'AGENCY]'));nnz(Agency_1)
Control_1 = NumBSC_1(contains(TextBSC_Env1(:,2),'CONTROL]'));nnz(Control_1)
Control2_1 = NumBSC_1(contains(TextBSC_Env1(:,2),'CONTROL2]'));nnz(Control2_1)

Threat_1 = NumBSC_1(contains(TextBSC_Env1(:,2),'THREAT]'));nnz(Threat_1)
Ownership_2 = NumBSC_2(contains(TextBSC_Env2(:,2),'OWNERSHIP]'));nnz(Ownership_2)
Agency_2 = NumBSC_2(contains(TextBSC_Env2(:,2),'AGENCY]'));nnz(Agency_2)
Control_2 = NumBSC_2(contains(TextBSC_Env2(:,2),'CONTROL]'));nnz(Control_2)
Control2_2 = NumBSC_2(contains(TextBSC_Env2(:,2),'CONTROL2]'));nnz(Control2_2)

Threat_2 = NumBSC_2(contains(TextBSC_Env2(:,2),'THREAT]'));nnz(Threat_2)


Ownership_Env1 = nanmean(Ownership_1); Agency_Env1 = nanmean(Agency_1); Control_Env1 = nanmean(Control_1); Control2_Env1 = nanmean(Control2_1);Threat_Env1 = nanmean(Threat_1);
Ownership_Env2 = nanmean(Ownership_2); Agency_Env2 = nanmean(Agency_2); Control_Env2 = nanmean(Control_2);Control2_Env2 = nanmean(Control2_2); Threat_Env2 = nanmean(Threat_2);

% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{'ParticipantsName'},{'OWNERSHIP_Env1'},{'AGENCY_Env1'},{'CONTROL_Env1'},{'CONTROL2_Env1'},{'THREAT_Env1'},{'OWNERSHIP_Env2'},{'AGENCY_Env2'},{'CONTROL_Env2'},{'CONTROL2_Env2'},{'THREAT_Env2'}],'BSC','A1:K1')
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',{MyParticipants(i,:)},'BSC',sprintf('A%s',num2str(i+1)))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[Ownership_Env1 Agency_Env1 Control_Env1 Control2_Env1 Threat_Env1 Ownership_Env2 Agency_Env2 Control_Env2 Control2_Env2 Threat_Env2],'BSC',sprintf('B%s:K%s',num2str(i+1),num2str(i+1)))
%  
%Timeline answers
if Correct_T_Env1(1) ~= 0
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{MyParticipants(i,:)}],'TimelineHit',sprintf('A%s',num2str(i+1)))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[T_Env1(Correct_T_Env1)'],'TimelineHit',sprintf('B%s:AH%s',num2str(i+1),num2str(i+1)))
end
if  Correct_T_Env2(1) ~= 0
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{MyParticipants(i,:)}],'TimelineHit',sprintf('A%s',num2str(i+1+size(MyParticipants,1))))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[T_Env2(Correct_T_Env2)'],'TimelineHit',sprintf('B%s:AH%s',num2str(i+1+size(MyParticipants,1)),num2str(i+1+size(MyParticipants,1))))
end
if Correct_F_Env1(1) ~= 0
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{MyParticipants(i,:)}],'TimelineCorrectRej',sprintf('A%s',num2str(i+1)))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[F_Env1(Correct_F_Env1)'],'TimelineCorrectRej',sprintf('B%s:AH%s',num2str(i+1),num2str(i+1)))
end
if  Correct_F_Env2(1) ~= 0
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{MyParticipants(i,:)}],'TimelineCorrectRej',sprintf('A%s',num2str(i+1+size(MyParticipants,1))))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[F_Env2(Correct_F_Env2)'],'TimelineCorrectRej',sprintf('B%s:AH%s',num2str(i+1+size(MyParticipants,1)),num2str(i+1+size(MyParticipants,1))))
end
if  Wrong_T_Env1(1)~= 0
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{MyParticipants(i,:)}],'TimelineMiss',sprintf('A%s',num2str(i+1)))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[T_Env1(Wrong_T_Env1)'],'TimelineMiss',sprintf('B%s:AH%s',num2str(i+1),num2str(i+1)))
end
if  Wrong_T_Env2(1) ~= 0
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{MyParticipants(i,:)}],'TimelineMiss',sprintf('A%s',num2str(i+1+size(MyParticipants,1))))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[T_Env2(Wrong_T_Env2)'],'TimelineMiss',sprintf('B%s:AH%s',num2str(i+1+size(MyParticipants,1)),num2str(i+1+size(MyParticipants,1))))
end
if Wrong_F_Env1~=0
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{MyParticipants(i,:)}],'TimelineFA',sprintf('A%s',num2str(i+1)))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[F_Env1(Wrong_F_Env1)'],'TimelineFA',sprintf('B%s:AH%s',num2str(i+1),num2str(i+1)))
end
if Wrong_F_Env2~=0
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{MyParticipants(i,:)}],'TimelineFA',sprintf('A%s',num2str(i+1+size(MyParticipants,1))))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[F_Env2(Wrong_F_Env2)'],'TimelineFA',sprintf('B%s:AH%s',num2str(i+1+size(MyParticipants,1)),num2str(i+1+size(MyParticipants,1))))
end
% figure(1)
% 
%     subplot(3,3,i)
%     if ConditionsAsynch(i) == 1 %if Env1 is asynch plot in black otherwise plot in blue
%         colors = 'ok'
%         colors2 = 'xb'
%     else
%         colors = 'ob'
%         colors2 = 'xk'
%     end
%     plot(Correct_T_Env1,ones(length(T_Env1(Correct_T_Env1))),colors)
%     hold on
%     plot(Correct_F_Env1,ones(length(F_Env1(Correct_F_Env1))),colors)
%     hold on
%        plot(Wrong_T_Env1,zeros(length(T_Env1(Wrong_T_Env1))),colors)
%     hold on
%     plot((Wrong_F_Env1),zeros(length(F_Env1(Wrong_F_Env1))),colors)
%     hold on    
%        plot(Correct_T_Env2,ones(length(T_Env2(Correct_T_Env2))),colors2)
%     hold on
%     plot(Correct_F_Env2,ones(length(F_Env2(Correct_F_Env2))),colors2)
%     hold on
%         if Wrong_T_Env2 ~= 0
%        plot((Wrong_T_Env2),zeros(length(T_Env2(Wrong_T_Env2))),colors2)
%         end
%     hold on
%     plot(Wrong_F_Env2,zeros(length(F_Env2(Wrong_F_Env2))),colors2)

% ylim([0 1.1])
%object difficulty
[F1_r1,F2_r1,F3_r1,F4_r1,F5_r1,F6_r1,F7_r1,F8_r1] = MistakeDetailed(sprintf('DesktopVRComputeranalysis/%s/%s',MyParticipants(i,:),ls(sprintf('DesktopVRComputeranalysis/%s/ROOMQ*.xlsx',MyParticipants(i,:)))),MyParticipants,i,'Env1') % check object difficulty for kitchen (env1)
[F1_r2,F2_r2,F3_r2,F4_r2,F5_r2,F6_r2,F7_r2,F8_r2] = MistakeDetailed(sprintf('DesktopVRComputeranalysis/%s/%s',MyParticipants(i,:),ls(sprintf('DesktopVRComputeranalysis/%s/ROOMQ*.xlsx',MyParticipants(i,:)))),MyParticipants,i,'Env2')% check object difficulty for labroom (env2)
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{MyParticipants(i,:)}],'Difficulty',sprintf('A%s',num2str(i+1)))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{'ParticipantsName'},{'F1_Env1'},{'F2_Env1'},{'F3_Env1'},{'F4_Env1'},{'F5_Env1'},{'F6_Env1'},{'F7_Env1'},{'F8_Env1'},{'F1_Env2'},{'F2_Env2'},{'F3_Env2'},{'F4_Env2'},{'F5_Env2'},{'F6_Env2'},{'F7_Env2'},{'F8_Env2'},{'zerofault_Env1'},{'zerofault_Env2'}],'Difficulty','A1:S1')
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[F1_r1 F2_r1 F3_r1 F4_r1 F5_r1 F6_r1 F7_r1 F8_r1 F1_r2 F2_r2 F3_r2 F4_r2 F5_r2 F6_r2 F7_r2 F8_r2 nnz([F1_r1 F2_r1 F3_r1 F4_r1 F5_r1 F6_r1 F7_r1 F8_r1]),nnz([F1_r2 F2_r2 F3_r2 F4_r2 F5_r2 F6_r2 F7_r2 F8_r2])],'Difficulty',sprintf('B%s:S%s',num2str(i+1),num2str(i+1)))
% 

% % Now lets have an excel sheet where each participant has each answer, and
% % then we check for each answer if this is correct or not.
% 
% %1. we do a vector of index containing both False and Catch trials
% Index_FalseAndCatch_ENV1 = IndexF_Env1;
% Index_FalseAndCatch_ENV1(find(Index_Catch_ENV1)) = 1;
% %2. We compute for Env1 and Env2 , whether the answer is correct or not
% for x = 1:length(NumQ1Env1)
%     if IndexT_Env1(x) == NumQ1Env1(x,1)
%         AnswerCorrectEnv1(x) = 1;
%      
%         else
%         
%         AnswerCorrectEnv1(x) = 0;
%  
%     end
% end
% CorrectRejEnv1 = NumQ1Env1(find(IndexT_Env1 == 1),1)
% CorrectRejEnv2 = NumQ1Env2(find(IndexT_Env2 == 1),1)
% HitEnv1 = NumQ1Env1(find(IndexT_Env1 == 0),1)
% HitEnv2 = NumQ1Env2(find(IndexT_Env2 == 0),1)
% 
% for x = 1:length(CorrectRejEnv1)
%     if IndexT_Env1(x) == CorrectRejEnv1(x)
%           AnswerEnv1CR(x) = 1;
%     else
%           AnswerEnv1CR(x) = 0;
%     end
% end
% 
% for x = 1:length(CorrectRejEnv2)
%     if IndexT_Env2(x) == CorrectRejEnv2(x)
%           AnswerEnv2CR(x) = 1;
%     else
%           AnswerEnv2CR(x) = 0;
%     end
% end
% for x = 1:length(HitEnv1)
%     if IndexT_Env1(x) == HitEnv1(x)
%           AnswerEnv1Hit(x) = 1;
%     else
%           AnswerEnv1Hit(x) = 0;
%     end
% end
% for x = 1:length(HitEnv2)
%     if IndexT_Env2(x) == HitEnv2(x)
%           AnswerEnv2Hit(x) = 1;
%     else
%           AnswerEnv2Hit(x) = 0;
%     end
% end
% for x = 1:length(NumQ1Env2)
%     if IndexT_Env2(i) == NumQ1Env2(x,1)
%         AnswerCorrectEnv2(x) = 1;
%        
%     else
%         AnswerCorrectEnv2(x) = 0;
%          
%     end
% end
% c = 1;
% % for z = 1:5:length(NumQ1Env1)
% %     if z+4<= length(NumQ1Env1)
% %     countEnv1(i,c) = nnz(AnswerCorrectEnv1(z:z+4));
% %     else
% %     countEnv1(i,c) = nnz(AnswerCorrectEnv1(z:end));
% %     end
% %     c =c+1;
% % end
% % c = 1;
% % for z = 1:5:length(NumQ1Env2)
% %    if z+4<= length(NumQ1Env2)
% %     countEnv2(i,c) = nnz(AnswerCorrectEnv2(z:z+4));
% %    else
% %     countEnv2(i,c) = nnz(AnswerCorrectEnv2(z:end));
% %    end
% %     c =c+1;
% % end
% c = 1
% for z = 1:length(NumQ1Env1)-4
%     if z+4<= length(NumQ1Env1)
%     countEnv1Window(i,c) = nnz(AnswerCorrectEnv1(z:z+4));
%     else
%     countEnv1Window(i,c) = nnz(AnswerCorrectEnv1(z:end));
%     end
%     c =c+1;
% end
% 
% c = 1
% for z = 1:length(CorrectRejEnv1)-4
%     if z+4<= length(CorrectRejEnv1)
%     countEnv1CRWindow(i,c) = nnz(AnswerEnv1CR(z:z+4));
%     else
%     countEnv1CRWindow(i,c) = nnz(AnswerEnv1CR(z:end));
%     end
%     c =c+1;
% end
% 
% c = 1
% for z = 1:length(CorrectRejEnv2)-4
%     if z+4<= length(CorrectRejEnv2)
%     countEnv2CRWindow(i,c) = nnz(AnswerEnv2CR(z:z+4));
%     else
%     countEnv2CRWindow(i,c) = nnz(AnswerEnv2CR(z:end));
%     end
%     c =c+1;
% end
% c = 1;
% for z = 1:length(HitEnv1)-4
%    if z+4<= length(HitEnv1)
%     countEnv1HitWindow(i,c) = nnz(AnswerEnv1Hit(z:z+4));
%    else
%     countEnv1HitWindow(i,c) = nnz(AnswerEnv1Hit(z:z+4));
%    end
%     c =c+1;
% end
% 
% c = 1;
% for z = 1:length(HitEnv2)-4
%    if z+4<= length(HitEnv2)
%     countEnv2HitWindow(i,c) = nnz(AnswerEnv2Hit(z:z+4));
%    else
%     countEnv2HitWindow(i,c) = nnz(AnswerEnv2Hit(z:z+4));
%    end
%     c =c+1;
% end
% 
% if list_env1issync(i) == 1 
%     countSynchWindow(i,1:46)= countEnv1Window(i,1:46);
%     countAsynchWindow(i,1:46) = countEnv2Window(i,1:46);
%     countSynchCRWindow(i,1:16)= countEnv1CRWindow(i,1:16);
%     countAsynchCRWindow(i,1:16) = countEnv2CRWindow(i,1:16);
%     countSynchHitWindow(i,1:26)= countEnv1HitWindow(i,1:26);
%     countAsynchHitWindow(i,1:26) = countEnv2HitWindow(i,1:26);
% else
%     countSynchWindow(i,1:46)= countEnv2Window(i,1:46);
%     countAsynchWindow(i,1:46) = countEnv1Window(i,1:46);
%     countSynchCRWindow(i,1:16)= countEnv2CRWindow(i,1:16);
%     countAsynchCRWindow(i,1:16) = countEnv1CRWindow(i,1:16);
%     countSynchHitWindow(i,1:26)= countEnv2HitWindow(i,1:26);
%     countAsynchHitWindow(i,1:26) = countEnv1HitWindow(i,1:26);
% end
% figure,
% plot(1:length( countSynchWindow(i,1:46)), countSynchWindow(i,1:46),'o')
% hold on
% if i == 22
% %we remove participant 10 and 30
%  countSynchWindow_R = countSynchWindow
%  countSynchWindow_R(1,:) = [];
%  countSynchWindow_R(12,:) = [];
%  countSynchCRWindow_R = countSynchCRWindow
%  countSynchCRWindow_R(1,:) = [];
%  countSynchCRWindow_R(12,:) = [];
%  countSynchHitWindow_R = countSynchHitWindow
%  countSynchHitWindow_R(1,:) = [];
%  countSynchHitWindow_R(12,:) = [];
%  countAsynchWindow_R = countAsynchWindow
%  countAsynchWindow_R(1,:) = [];
%  countAsynchWindow_R(12,:) = [];
%  countAsynchCRWindow_R = countAsynchCRWindow
%  countAsynchCRWindow_R(1,:) = [];
%  countAsynchCRWindow_R(12,:) = [];
%  countAsynchHitWindow_R = countAsynchHitWindow
%  countAsynchHitWindow_R(1,:) = [];
%  countAsynchHitWindow_R(12,:) = [];
% %now in term of cumulative 
% for u = 1:size(MyParticipants,1)-2
% for l = 1:size(countSynchWindow,2)
% CumulativeCountSynch_RWindow(u,l) = sum(countSynchWindow_R(u,1:l))
% end
% for b = 1:size(countSynchCRWindow_R,2)
% CumulativeCountSynchCR_RWindow(u,b) = sum(countSynchCRWindow_R(u,1:b))
% CumulativeCountAsynchCR_RWindow(u,b) = sum(countAsynchCRWindow_R(u,1:b))
% end
% for v = 1:size(countSynchHitWindow_R,2)
% CumulativeCountSynchHit_RWindow(u,v) = sum(countSynchHitWindow_R(u,1:v))
% CumulativeCountAsynchHit_RWindow(u,v) = sum(countAsynchHitWindow_R(u,1:v))
% end
% end
% end
% if i ==22
%   for k = 1:length(MyParticipants)-2
% 
%   figure(9)
%     subplot(1,2,1)
%     title('Synch')
%     ylabel('cumulative correct answer over 5 trials')
%     plot(1:size(CumulativeCountSynch_RWindow,2),CumulativeCountSynch_RWindow(k,:),'-');
%     hold on;
%     if k == length(MyParticipants)-2
%     plot(1:size(CumulativeCountSynch_RWindow,2),mean(CumulativeCountSynch_RWindow,1),'-k','LineWidth',2)
%     hold on;
%     end
%     subplot(1,2,2)
%     title('Asynch')
%     plot(1:size(CumulativeCountAsynch_RWindow,2),CumulativeCountAsynch_RWindow(k,:),'-');
%     hold on;
%     if k == length(MyParticipants)-2
%     plot(1:size(CumulativeCountAsynch_RWindow,2),mean(CumulativeCountAsynch_RWindow,1),'-k','LineWidth',2)
%     hold on;
%     end
%   end
% 
% 
% %CR
%   for k = 1:length(MyParticipants)-2
% 
%   figure(18)
%     subplot(1,2,1)
%     title('CR Synch')
%     ylabel('cumulative correct rejection')
%     plot(1:size(CumulativeCountSynchCR_RWindow,2),CumulativeCountSynchCR_RWindow(k,:),'-');
%     hold on;
%     if k == length(MyParticipants)-2
%     plot(1:size(CumulativeCountSynchCR_RWindow,2),mean(CumulativeCountSynchCR_RWindow,1),'-k','LineWidth',2)
%     hold on;
%     end
%     subplot(1,2,2)
%     title('CR Asynch')
%     plot(1:size(CumulativeCountAsynchCR_RWindow,2),CumulativeCountAsynchCR_RWindow(k,:),'-');
%     hold on;
%     if k == length(MyParticipants)-2
%     plot(1:size(CumulativeCountAsynchCR_RWindow,2),mean(CumulativeCountAsynchCR_RWindow,1),'-k','LineWidth',2)
%     hold on;
%     end
%   end
%   
%   %Hit
%   for k = 1:length(MyParticipants)-2
% 
%   figure(19)
%     subplot(1,2,1)
%     title('Hit Synch')
%     ylabel('cumulative correct hit')
%     plot(1:size(CumulativeCountSynchHit_RWindow,2),CumulativeCountSynchHit_RWindow(k,:),'-');
%     hold on;
%     if k == length(MyParticipants)-2
%     plot(1:size(CumulativeCountSynchHit_RWindow,2),mean(CumulativeCountSynchHit_RWindow,1),'-k','LineWidth',2)
%     hold on;
%     end
%     subplot(1,2,2)
%     title('Hit Asynch')
%     plot(1:size(CumulativeCountAsynchHit_RWindow,2),CumulativeCountAsynchHit_RWindow(k,:),'-');
%     hold on;
%     if k == length(MyParticipants)-2
%     plot(1:size(CumulativeCountAsynchHit_RWindow,2),mean(CumulativeCountAsynchHit_RWindow,1),'-k','LineWidth',2)
%     hold on;
%     end
%   end
% % now we fit a line per participants, we extract the slope
% for u = 1:size(MyParticipants,1)-2
% coefficients = polyfit((1:46), CumulativeCountSynch_RWindow(u,:), 1);
% SlopeSynch(u) = coefficients(1);
% coefficientsA = polyfit((1:46), CumulativeCountAsynch_RWindow(u,:), 1);
% SlopeAsynch(u) = coefficientsA(1);
% 
% coefficientsCR = polyfit((1:16), CumulativeCountSynchCR_RWindow(u,:), 1);
% SlopeSynchCR(u) = coefficientsCR(1);
% coefficientsACR = polyfit((1:16), CumulativeCountAsynchCR_RWindow(u,:), 1);
% SlopeAsynchCR(u) = coefficientsACR(1);
% 
% coefficientsHit = polyfit((1:26), CumulativeCountSynchHit_RWindow(u,:), 1);
% SlopeSynchHit(u) = coefficientsHit(1);
% coefficientsAHit = polyfit((1:26), CumulativeCountAsynchHit_RWindow(u,:), 1);
% SlopeAsynchHit(u) = coefficientsAHit(1);
% end
% [h,p] = ttest(SlopeSynch,SlopeAsynch)
% [hcr,pcr] = ttest(SlopeSynchCR,SlopeAsynchCR)
% [hhit,phit] = ttest(SlopeSynchHit,SlopeAsynchHit)
% end
% g = 1;
% for u = 1:size(MyParticipants)-2
%     for j = 1:size(CumulativeCountSynch_RWindow,2)-1
% SlopeSynchPerPointPerParticipant(u,g) = CumulativeCountSynch_RWindow(u,j+1)-CumulativeCountSynch_RWindow(u,j);
% SlopeAsynchPerPointPerParticipant(u,g) = CumulativeCountAsynch_RWindow(u,j+1)-CumulativeCountAsynch_RWindow(u,j);
% g = g+1;
%     end
% end

%now the same but per step of 10 instead

%  for u = 1:size(MyParticipants,1)
%      t = 1;
% for z = 1:2:10
%    
%     countSynch_10(u,t)= countSynch(u,z)+countSynch(u,z+1);
%     countAsynch_10(u,t) = countAsynch(u,z)+countAsynch(u,z+1);
%     t=t+1;
%     end
%  end
% 
%  
 %now in term of cumulative 
% for u = 1:size(MyParticipants,1)-2
% for l = 1:5
% CumulativeCountSynch_R(u,l) = sum(countSynch_R(u,1:l))
% CumulativeCountAsynch_R(u,l) = sum(countAsynch_R(u,1:l))
% end
% end
% 
% for u = 1:size(MyParticipants,1)-2
% for l = 1:10
% CumulativeCountSynch_R5(u,l) = sum(countSynch_R5(u,1:l))
% CumulativeCountAsynch_R5(u,l) = sum(countAsynch_R5(u,1:l))
% end
% end
% if k == length(MyParticipants)
%     save countSynch;
%     save countAsynch;
% end
% %we remove participant 10 and 30
%  countSynch_R5 = countSynch
%  countSynch_R5(1,:) = [];
%  countSynch_R5(12,:) = [];
%  
%  countAsynch_R5 = countAsynch
%  countAsynch_R5(1,:) = [];
%  countAsynch_R5(12,:) = [];
% 
% if i == length(MyParticipants)
%     for k = 1:length(MyParticipants)-2
%     figure(1)
%     subplot(1,2,1)
%     title('Synch')
%     ylabel('correct answer over 5 trials')
%     plot(1:size(countSynch_R5,2),countSynch_R5(k,:),'-o');
%     hold on;
%     if k == length(MyParticipants)-2
%     plot(1:size(countSynch_R5,2),mean(countSynch_R5,1),'-ok','LineWidth',2)
%     hold on;
%     end
%     subplot(1,2,2)
%     title('Asynch')
%     plot(1:size(countAsynch_R5,2),countAsynch_R5(k,:),'-o');
%     hold on;
%     if k == length(MyParticipants)-2
%     plot(1:size(countAsynch_R5,2),mean(countAsynch_R5,1),'-ok','LineWidth',2)
%     hold on;
%     end
%     end
% end
% 
% %we remove participant 10 and 30
%  countSynch_R = countSynch_10
%  countSynch_R(1,:) = [];
%  countSynch_R(12,:) = [];
%  
%  countAsynch_R = countAsynch_10
%  countAsynch_R(1,:) = [];
%  countAsynch_R(12,:) = [];
% if i == length(MyParticipants)
%     for k = 1:length(MyParticipants)-2
%     figure(2)
%     subplot(1,2,1)
%     title('Synch')
%     plot(1:size(countSynch_R,2),countSynch_R(k,:),'-o');
%     hold on;
%     ylabel('correct answer over 10 trials')
%     if k == length(MyParticipants)-2
%     plot(1:size(countSynch_10,2),mean(countSynch_R,1),'-ok','LineWidth',2)
%     hold on;
%     end
%     subplot(1,2,2)
%     plot(1:size(countAsynch_R,2),countAsynch_R(k,:),'-o');
%     hold on;
%     if k == length(MyParticipants)-2
%     title('Asynch')
%     plot(1:size(countAsynch_R,2),mean(countAsynch_R,1),'-ok','LineWidth',2)
%     hold on;
%     end
%     end
% end
% 
% 
% %with cumulative rate
% %per five trials
% 
%  for k = 1:length(MyParticipants)-2
%     figure(5)
%     subplot(1,2,1)
%     title('Synch')
%     ylabel('cumulative correct answer over 5 trials')
%     plot(1:size(CumulativeCountSynch_R5,2),CumulativeCountSynch_R5(k,:),'-o');
%     hold on;
%     if k == length(MyParticipants)-2
%     plot(1:size(CumulativeCountSynch_R5,2),mean(CumulativeCountSynch_R5,1),'-ok','LineWidth',2)
%     hold on;
%     end
%     subplot(1,2,2)
%     title('Asynch')
%     plot(1:size(CumulativeCountAsynch_R5,2),CumulativeCountAsynch_R5(k,:),'-o');
%     hold on;
%     if k == length(MyParticipants)-2
%     plot(1:size(CumulativeCountAsynch_R5,2),mean(CumulativeCountAsynch_R5,1),'-ok','LineWidth',2)
%     hold on;
%     end
%  end
%     
%  
%  %per 10 trials
%   for k = 1:length(MyParticipants)-2
% 
%   figure(6)
%     subplot(1,2,1)
%     title('Synch')
%     ylabel('cumulative correct answer over 5 trials')
%     plot(1:size(CumulativeCountSynch_R,2),CumulativeCountSynch_R(k,:),'-o');
%     hold on;
%     if k == length(MyParticipants)-2
%     plot(1:size(CumulativeCountSynch_R,2),mean(CumulativeCountSynch_R,1),'-ok','LineWidth',2)
%     hold on;
%     end
%     subplot(1,2,2)
%     title('Asynch')
%     plot(1:size(CumulativeCountAsynch_R,2),CumulativeCountAsynch_R(k,:),'-o');
%     hold on;
%     if k == length(MyParticipants)-2
%     plot(1:size(CumulativeCountAsynch_R,2),mean(CumulativeCountAsynch_R,1),'-ok','LineWidth',2)
%     hold on;
%     end
%   end
%3. on écrit dans un xlsfile
% if list_env1issync(i) == 1 
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{MyParticipants(i,:)}],'AnswerSynch',sprintf('A%s:A%s',num2str(i+1+ (length(50)*(i-1))),num2str(i+1+ (length(50)*(i-1))+length(NumQ1Env1))))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{'ParticipantsName'},{'Answer'},{'Correct'},{'CorrectAnswer'}],'AnswerSynch','A1:D1')
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[NumQ1Env1(:,1) IndexT_Env1 AnswerCorrectEnv1'],'AnswerSynch',sprintf('B%s:D%s',num2str(i+1+ (length(50)*(i-1))),num2str(i+1+ (length(50)*(i-1))+length(NumQ1Env1))))
% 
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{MyParticipants(i,:)}],'AnswerAsynch',sprintf('A%s:A%s',num2str(i+1+ (length(50)*(i-1))),num2str(i+1+ (length(50)*(i-1))+length(NumQ1Env2))))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{'ParticipantsName'},{'Answer'},{'Correct'},{'CorrectAnswer'}],'AnswerAsynch','A1:D1')
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[NumQ1Env2(:,1) IndexT_Env2 AnswerCorrectEnv2'],'AnswerAsynch',sprintf('B%s:D%s',num2str(i+1+ (length(50)*(i-1))),num2str(i+1+ (length(50)*(i-1))+length(NumQ1Env2))))
% 
% else if list_env1issync(i) == 2
%         %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{MyParticipants(i,:)}],'AnswerAsynch',sprintf('A%s:A%s',num2str(i+1+ (length(50)*(i-1))),num2str(i+1+ (length(50)*(i-1))+length(NumQ1Env1))))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{'ParticipantsName'},{'Answer'},{'Correct'},{'CorrectAnswer'}],'AnswerAsynch','A1:D1')
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[NumQ1Env1(:,1) IndexT_Env1 AnswerCorrectEnv1'],'AnswerAsynch',sprintf('B%s:D%s',num2str(i+1+ (length(50)*(i-1))),num2str(i+1+ (length(50)*(i-1))+length(NumQ1Env1))))
% 
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{MyParticipants(i,:)}],'AnswerSynch',sprintf('A%s:A%s',num2str(i+1+ (length(50)*(i-1))),num2str(i+1+ (length(50)*(i-1))+length(NumQ1Env2))))
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[{'ParticipantsName'},{'Answer'},{'Correct'},{'CorrectAnswer'}],'AnswerSynch','A1:D1')
% %xlswrite('Analysis_Participants_ABMPilotV2_ALLANSWERs.xlsx',[NumQ1Env2(:,1) IndexT_Env2 AnswerCorrectEnv2'],'AnswerSynch',sprintf('B%s:D%s',num2str(i+1+ (length(50)*(i-1))),num2str(i+1+ (length(50)*(i-1))+length(NumQ1Env2))))
% 
%     end
% end
% figure(2)
% if ConditionsAsynch(i) == 1 %if Env1 is asynch plot in black otherwise plot in blue
%         colors = 'ok'
%         colors2 = 'ob'
%     else
%         colors = 'ob'
%         colors2 = 'ok'
%     end
% subplot(4,3,i)
% plot(1:length(NumQ2(IndexEnv1)),NumQ2(IndexEnv1),colors)
% hold on
% plot(1:length(NumQ2(IndexEnv2)),NumQ2(IndexEnv2),colors2)
% hold on


 %Num = []; text = []; All = [];NumQ1 = []; NumQ2 = []; NumQ1Env1 = []; NumQ2Env1 = [];NumQ1Env2 = []; NumQ2Env2 = [];AnswerCorrectEnv1 = [];AnswerCorrectEnv2 = []
end