%% clear everything
clc; clear; close all

%%% Group/Condition one path of data folder and group name
path='/Users/snray/Documents/data_lab/AN/BoldData/';
cd(path)
SUBJlist_Group1=dir('AN*');
group ='pre';
%% Absolute and Random and Normalised CC, PL, SW for Group1 data extraction
for i = 1:length(SUBJlist_Group1)
    %%
    SUBJname=SUBJlist_Group1(i).name
    cd(SUBJname)
    cd(group)

    Group1fix_sub_name=SUBJname; %14
    data=load([Group1fix_sub_name '_pre_ABS.mat'])
    Group1_CC(i,:,:)= data.GT_clust_coeff;         
    Group1_PL(i,:,:)= data.GT_path_length;
    Group1_LE(i,:,:)= data.GT_local_eff;         
    Group1_GE(i,:,:)= data.GT_global_eff;
    Group1_Degree(i,:,:)= data.GT_degree;
    Group1_PC(i,:,:)= data.GT_participation_coeff;  
    Group1_Mod(i,:,:)= data.GT_modularity; 
    Group1_ComStr(i,:,:)= data.GT_community_structure; 
    Group1_BC(i,:,:)= data.GT_betw_cent; 
    Group1_EC(i,:,:)= data.GT_eigv_cent; 
    Group1_FC(i,:,:)= data.GT_flow_coef; 
    Group1_AS(i,:,:)= data.GT_assortvity;
    Group1_KC(i,:,:)= data.GT_k_core; 
    
    data_rand=load([Group1fix_sub_name '_pre_RAND.mat'])
    Group1_CC_rand(i,:,:,:)= data_rand.GT_clust_coeff_rand;     
    Group1_PL_rand(i,:,:,:)= data_rand.GT_path_length_rand;  
    Group1_LE_rand(i,:,:,:)= data_rand.GT_local_eff_rand;    
    Group1_GE_rand(i,:,:,:)= data_rand.GT_global_eff_rand ;  %GT_global_eff_rand;
    Group1_Degree_rand(i,:,:,:)= data_rand.GT_degree_rand;
    Group1_PC_rand(i,:,:,:)= data_rand.GT_participation_coeff_rand;
    Group1_Mod_rand(i,:,:,:)= data_rand.GT_modularity_rand; 
    Group1_ComStr_rand(i,:,:,:)= data_rand.GT_community_structure_rand; 
    Group1_BC_rand(i,:,:,:)= data_rand.GT_betw_cent_rand; 
    Group1_EC_rand(i,:,:,:)= data_rand.GT_eigv_cent_rand; 
    Group1_FC_rand(i,:,:,:)= data_rand.GT_flow_coef_rand; 
    Group1_AS_rand(i,:,:,:)= data_rand.GT_assortvity_rand;
    Group1_KC_rand(i,:,:,:)= data_rand.GT_k_core_rand;   
    cd ..
    cd ..
end
%%
SparcityStart=1;
SparcityEnd=length(data.sparsity_val);
NoOfRoi=size(Group1_CC,3);
%%
Group1_CC_50=Group1_CC(:,SparcityStart:SparcityEnd,:);                          
Group1_CC_rand_squ = squeeze(mean(Group1_CC_rand,3));        
AvgGroup1_CC=mean(mean(Group1_CC_50,3));

sparsity_CC_Group1_50 = (mean(Group1_CC_50,3));
sparsity_CC_rand_Group1 = mean(Group1_CC_rand_squ,3);
sparsity_CC_rand_Group1_50 = sparsity_CC_rand_Group1 (:,SparcityStart:SparcityEnd);

Group1_PL_50=Group1_PL(:,:,SparcityStart:SparcityEnd);
%AvgGroup1_PL=squeeze(mean(Group1_PL_50));
Sparsity_PL_Group1_50=squeeze(Group1_PL_50);
Sparsity_PL_Group1_rand_squ = squeeze(mean(Group1_PL_rand,3)); 
Sparsity_PL_Group1_rand_50 = Sparsity_PL_Group1_rand_squ (:,SparcityStart:SparcityEnd);

Group1_Degree_50=Group1_Degree(:,SparcityStart:SparcityEnd,:);                             
AvgGroup1_Degree=mean(mean(Group1_Degree_50,3));
sparsity_Degree_Group1_50 = (mean(Group1_Degree_50,3));

Group1_PC_50=Group1_PC(:,SparcityStart:SparcityEnd,:);                          
Group1_PC_rand_squ = squeeze(mean(Group1_PC_rand,3));        
AvgGroup1_PC=mean(mean(Group1_PC_50,3));

sparsity_PC_Group1_50 = (mean(Group1_PC_50,3));
sparsity_PC_rand_Group1 = mean(Group1_PC_rand_squ,3);
sparsity_PC_rand_Group1_50 = sparsity_PC_rand_Group1 (:,SparcityStart:SparcityEnd);

Group1_BC_50=Group1_BC(:,SparcityStart:SparcityEnd,:);                          
Group1_BC_rand_squ = squeeze(mean(Group1_BC_rand,3));        
sparsity_BC_Group1_50 = (mean(Group1_BC_50,3));
sparsity_BC_rand_Group1 = mean(Group1_BC_rand_squ,3);
sparsity_BC_rand_Group1_50 = sparsity_BC_rand_Group1 (:,SparcityStart:SparcityEnd);

Group1_EC_50=Group1_EC(:,SparcityStart:SparcityEnd,:);                          
Group1_EC_rand_squ = squeeze(mean(Group1_EC_rand,3));        
sparsity_EC_Group1_50 = (mean(Group1_EC_50,3));
sparsity_EC_rand_Group1 = mean(Group1_EC_rand_squ,3);
sparsity_EC_rand_Group1_50 = sparsity_EC_rand_Group1 (:,SparcityStart:SparcityEnd);

Group1_FC_50=Group1_FC(:,SparcityStart:SparcityEnd,:);                          
Group1_FC_rand_squ = squeeze(mean(Group1_FC_rand,3));        
sparsity_FC_Group1_50 = (mean(Group1_FC_50,3));
sparsity_FC_rand_Group1 = mean(Group1_FC_rand_squ,3);
sparsity_FC_rand_Group1_50 = sparsity_FC_rand_Group1 (:,SparcityStart:SparcityEnd);

Group1_KC_50=Group1_KC(:,SparcityStart:SparcityEnd,:);                          
Group1_KC_rand_squ = squeeze(mean(Group1_KC_rand,3));        
sparsity_KC_Group1_50 = (mean(Group1_KC_50,3));
sparsity_KC_rand_Group1 = mean(Group1_KC_rand_squ,3);
sparsity_KC_rand_Group1_50 = sparsity_KC_rand_Group1 (:,SparcityStart:SparcityEnd);

Group1_Mod_50=Group1_Mod(:,:,SparcityStart:SparcityEnd);
Sparsity_Mod_Group1_50=squeeze(Group1_Mod_50);
Sparsity_Mod_Group1_rand_squ = squeeze(mean(Group1_Mod_rand,3)); 
Sparsity_Mod_Group1_rand_50 = Sparsity_Mod_Group1_rand_squ (:,SparcityStart:SparcityEnd);

Group1_AS_50=Group1_AS(:,:,SparcityStart:SparcityEnd);
Sparsity_AS_Group1_50=squeeze(Group1_AS_50);
Sparsity_AS_Group1_rand_squ = squeeze(mean(Group1_AS_rand,3)); 
Sparsity_AS_Group1_rand_50 = Sparsity_AS_Group1_rand_squ (:,SparcityStart:SparcityEnd);
%%
for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_CC_normalised_Group1(i,j) = sparsity_CC_Group1_50(i,j)/sparsity_CC_rand_Group1_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_PL_normalised_Group1(i,j) = Sparsity_PL_Group1_50(i,j)/Sparsity_PL_Group1_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            SmallWorldNess_Group1(i,j) = sparsity_CC_normalised_Group1(i,j)/sparsity_PL_normalised_Group1(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_PC_normalised_Group1(i,j) = sparsity_PC_Group1_50(i,j)/sparsity_PC_rand_Group1_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_BC_normalised_Group1(i,j) = sparsity_BC_Group1_50(i,j)/sparsity_BC_rand_Group1_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_EC_normalised_Group1(i,j) = sparsity_EC_Group1_50(i,j)/sparsity_EC_rand_Group1_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_FC_normalised_Group1(i,j) = sparsity_FC_Group1_50(i,j)/sparsity_FC_rand_Group1_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_KC_normalised_Group1(i,j) = sparsity_KC_Group1_50(i,j)/sparsity_KC_rand_Group1_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_Mod_normalised_Group1(i,j) = Sparsity_Mod_Group1_50(i,j)/Sparsity_Mod_Group1_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_AS_normalised_Group1(i,j) = Sparsity_AS_Group1_50(i,j)/Sparsity_AS_Group1_rand_50(i,j); 
    end
end
%%
Group1_LE_50=Group1_LE(:,SparcityStart:SparcityEnd,:);                          
Group1_LE_rand_squ = squeeze(mean(Group1_LE_rand,3));        
AvgGroup1_LE=mean(mean(Group1_LE_50,3));

sparsity_LE_Group1_50 = (mean(Group1_LE_50,3));
sparsity_LE_rand_Group1 = mean(Group1_LE_rand_squ,3);
sparsity_LE_rand_Group1_50 = sparsity_LE_rand_Group1 (:,SparcityStart:SparcityEnd);

Group1_GE_50=Group1_GE(:,SparcityStart:SparcityEnd);
AvgGroup1_GE=squeeze(mean(Group1_GE_50));

Sparsity_GE_Group1_50=squeeze(Group1_GE_50);
Sparsity_GE_Group1_rand_squ = squeeze(mean(Group1_GE_rand,3)); 
Sparsity_GE_Group1_rand_50 = Sparsity_GE_Group1_rand_squ (:,SparcityStart:SparcityEnd); 

%%
for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_LE_normalised_Group1(i,j) = sparsity_LE_Group1_50(i,j)/sparsity_LE_rand_Group1_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_GE_normalised_Group1(i,j) = Sparsity_GE_Group1_50(i,j)/Sparsity_GE_Group1_rand_50(i,j); 
    end
end




%% Group/Condition two path of data folder and group name
SUBJlist_Group2=dir('AN*');
group ='post';
%%% Absolute and Random and Normalised CC, PL, SW for Group2 data extraction
for i = 1:length(SUBJlist_Group2)
    %%
    SUBJname=SUBJlist_Group2(i).name
    cd(SUBJname)
    cd(group)
    Group1fix_sub_name=SUBJname; %14
    data=load([Group1fix_sub_name '_post_ABS.mat'])
    Group2_CC(i,:,:)= data.GT_clust_coeff;         
    Group2_PL(i,:,:)= data.GT_path_length;
    Group2_LE(i,:,:)= data.GT_local_eff;         
    Group2_GE(i,:,:)= data.GT_global_eff;
    Group2_Degree(i,:,:)= data.GT_degree;
    Group2_PC(i,:,:)= data.GT_participation_coeff;
    Group2_Mod(i,:,:)= data.GT_modularity; 
    Group2_ComStr(i,:,:)= data.GT_community_structure; 
    Group2_BC(i,:,:)= data.GT_betw_cent; 
    Group2_EC(i,:,:)= data.GT_eigv_cent; 
    Group2_FC(i,:,:)= data.GT_flow_coef; 
    Group2_AS(i,:,:)= data.GT_assortvity;
    Group2_KC(i,:,:)= data.GT_k_core; 
    
    data_rand=load([Group1fix_sub_name '_post_RAND.mat'])
    Group2_CC_rand(i,:,:,:)= data_rand.GT_clust_coeff_rand;     
    Group2_PL_rand(i,:,:,:)= data_rand.GT_path_length_rand;  
    Group2_LE_rand(i,:,:,:)= data_rand.GT_local_eff_rand;    
    Group2_GE_rand(i,:,:,:)= data_rand.GT_global_eff_rand;
    Group2_Degree_rand(i,:,:,:)= data_rand.GT_degree_rand;
    Group2_PC_rand(i,:,:,:)= data_rand.GT_participation_coeff_rand;
    Group2_Mod_rand(i,:,:,:)= data_rand.GT_modularity_rand; 
    Group2_ComStr_rand(i,:,:,:)= data_rand.GT_community_structure_rand; 
    Group2_BC_rand(i,:,:,:)= data_rand.GT_betw_cent_rand; 
    Group2_EC_rand(i,:,:,:)= data_rand.GT_eigv_cent_rand; 
    Group2_FC_rand(i,:,:,:)= data_rand.GT_flow_coef_rand; 
    Group2_AS_rand(i,:,:,:)= data_rand.GT_assortvity_rand;
    Group2_KC_rand(i,:,:,:)= data_rand.GT_k_core_rand;  
    cd ..
    cd ..
end
%%
Group2_CC_50=Group2_CC(:,SparcityStart:SparcityEnd,:);                          
Group2_CC_rand_squ = squeeze(mean(Group2_CC_rand,3));        
AvgGroup2_CC=mean(mean(Group2_CC_50,3));

sparsity_CC_Group2_50 = (mean(Group2_CC_50,3));
sparsity_CC_rand_Group2 = mean(Group2_CC_rand_squ,3);
sparsity_CC_rand_Group2_50 = sparsity_CC_rand_Group2 (:,SparcityStart:SparcityEnd);

Group2_PL_50=Group2_PL(:,:,SparcityStart:SparcityEnd);
%AvgGroup2_PL=squeeze(mean(Group2_PL_50));
Sparsity_PL_Group2_50=squeeze(Group2_PL_50);
Sparsity_PL_Group2_rand_squ = squeeze(mean(Group2_PL_rand,3)); 
Sparsity_PL_Group2_rand_50 = Sparsity_PL_Group2_rand_squ (:,SparcityStart:SparcityEnd);

Group2_Degree_50=Group2_Degree(:,SparcityStart:SparcityEnd,:);                             
AvgGroup2_Degree=mean(mean(Group2_Degree_50,3));
sparsity_Degree_Group2_50 = (mean(Group2_Degree_50,3));

Group2_PC_50=Group2_PC(:,SparcityStart:SparcityEnd,:);                          
Group2_PC_rand_squ = squeeze(mean(Group2_PC_rand,3));        
AvgGroup2_PC=mean(mean(Group2_PC_50,3));

sparsity_PC_Group2_50 = (mean(Group2_PC_50,3));
sparsity_PC_rand_Group2 = mean(Group2_PC_rand_squ,3);
sparsity_PC_rand_Group2_50 = sparsity_PC_rand_Group2 (:,SparcityStart:SparcityEnd);

Group2_BC_50=Group2_BC(:,SparcityStart:SparcityEnd,:);                          
Group2_BC_rand_squ = squeeze(mean(Group2_BC_rand,3));        
sparsity_BC_Group2_50 = (mean(Group2_BC_50,3));
sparsity_BC_rand_Group2 = mean(Group2_BC_rand_squ,3);
sparsity_BC_rand_Group2_50 = sparsity_BC_rand_Group2 (:,SparcityStart:SparcityEnd);

Group2_EC_50=Group2_EC(:,SparcityStart:SparcityEnd,:);                          
Group2_EC_rand_squ = squeeze(mean(Group2_EC_rand,3));        
sparsity_EC_Group2_50 = (mean(Group2_EC_50,3));
sparsity_EC_rand_Group2 = mean(Group2_EC_rand_squ,3);
sparsity_EC_rand_Group2_50 = sparsity_EC_rand_Group2 (:,SparcityStart:SparcityEnd);

Group2_FC_50=Group2_FC(:,SparcityStart:SparcityEnd,:);                          
Group2_FC_rand_squ = squeeze(mean(Group2_FC_rand,3));        
sparsity_FC_Group2_50 = (mean(Group2_FC_50,3));
sparsity_FC_rand_Group2 = mean(Group2_FC_rand_squ,3);
sparsity_FC_rand_Group2_50 = sparsity_FC_rand_Group2 (:,SparcityStart:SparcityEnd);

Group2_KC_50=Group2_KC(:,SparcityStart:SparcityEnd,:);                          
Group2_KC_rand_squ = squeeze(mean(Group2_KC_rand,3));        
sparsity_KC_Group2_50 = (mean(Group2_KC_50,3));
sparsity_KC_rand_Group2 = mean(Group2_KC_rand_squ,3);
sparsity_KC_rand_Group2_50 = sparsity_KC_rand_Group2 (:,SparcityStart:SparcityEnd);

Group2_Mod_50=Group2_Mod(:,:,SparcityStart:SparcityEnd);
Sparsity_Mod_Group2_50=squeeze(Group2_Mod_50);
Sparsity_Mod_Group2_rand_squ = squeeze(mean(Group2_Mod_rand,3)); 
Sparsity_Mod_Group2_rand_50 = Sparsity_Mod_Group2_rand_squ (:,SparcityStart:SparcityEnd);

Group2_AS_50=Group2_AS(:,:,SparcityStart:SparcityEnd);
Sparsity_AS_Group2_50=squeeze(Group2_AS_50);
Sparsity_AS_Group2_rand_squ = squeeze(mean(Group2_AS_rand,3)); 
Sparsity_AS_Group2_rand_50 = Sparsity_AS_Group2_rand_squ (:,SparcityStart:SparcityEnd);
%%
for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_CC_normalised_Group2(i,j) = sparsity_CC_Group2_50(i,j)/sparsity_CC_rand_Group2_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_PL_normalised_Group2(i,j) = Sparsity_PL_Group2_50(i,j)/Sparsity_PL_Group2_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            SmallWorldNess_Group2(i,j) = sparsity_CC_normalised_Group2(i,j)/sparsity_PL_normalised_Group2(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_PC_normalised_Group2(i,j) = sparsity_PC_Group2_50(i,j)/sparsity_PC_rand_Group2_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_BC_normalised_Group2(i,j) = sparsity_BC_Group2_50(i,j)/sparsity_BC_rand_Group2_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_EC_normalised_Group2(i,j) = sparsity_EC_Group2_50(i,j)/sparsity_EC_rand_Group2_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_FC_normalised_Group2(i,j) = sparsity_FC_Group2_50(i,j)/sparsity_FC_rand_Group2_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_KC_normalised_Group2(i,j) = sparsity_KC_Group2_50(i,j)/sparsity_KC_rand_Group2_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_Mod_normalised_Group2(i,j) = Sparsity_Mod_Group2_50(i,j)/Sparsity_Mod_Group2_rand_50(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_AS_normalised_Group2(i,j) = Sparsity_AS_Group2_50(i,j)/Sparsity_AS_Group2_rand_50(i,j); 
    end
end
%%
Group2_LE_50=Group2_LE(:,SparcityStart:SparcityEnd,:);                          
Group2_LE_rand_squ = squeeze(mean(Group2_LE_rand,3));        
AvgGroup2_LE=mean(mean(Group2_LE_50,3));

sparsity_LE_Group2_50 = (mean(Group2_LE_50,3));
sparsity_LE_rand_Group2 = mean(Group2_LE_rand_squ,3);
sparsity_LE_rand_Group2_50 = sparsity_LE_rand_Group2 (:,SparcityStart:SparcityEnd);

Group2_GE_50=Group2_GE(:,SparcityStart:SparcityEnd);
AvgGroup2_GE=squeeze(mean(Group2_GE_50));

Sparsity_GE_Group2_50=squeeze(Group2_GE_50);
Sparsity_GE_Group2_rand_squ = squeeze(mean(Group2_GE_rand,3)); 
Sparsity_GE_Group2_rand_50 = Sparsity_GE_Group2_rand_squ (:,SparcityStart:SparcityEnd); 

%%
for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_LE_normalised_Group2(i,j) = sparsity_LE_Group2_50(i,j)/sparsity_LE_rand_Group2_50(i,j); 
    end
end


for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_GE_normalised_Group2(i,j) = Sparsity_GE_Group2_50(i,j)/Sparsity_GE_Group2_rand_50(i,j); 
    end
end

%% Absolute CC and PC ploting
ci=1.5; % 1.64; %99.5=2.807, 99=2.576, 95=1.960, 90=1.645, 85=1.440, 80=1.282, ci=1 is SEM (Standard Error Of The Mean)
figure;
y1_CC = mean (sparsity_CC_Group1_50);
z1_CC = std (sparsity_CC_Group1_50)/sqrt (length (sparsity_CC_Group1_50)); 
errorbar (y1_CC,ci*z1_CC, 'r','LineWidth',1); grid on; 
hold on
y2_CC = mean (sparsity_CC_Group2_50);
z2_CC = std (sparsity_CC_Group2_50)/sqrt (length (sparsity_CC_Group2_50));
errorbar (y2_CC,ci*z2_CC, 'b','LineWidth',1); grid on; 
legend('Pre TMS', 'Post TMS')
title('Absolute Clustering Coefficient')
figure;
y1_PC = mean (sparsity_PC_Group1_50);
z1_PC = std (sparsity_PC_Group1_50)/sqrt (length (sparsity_PC_Group1_50)); 
errorbar (y1_PC,ci*z1_PC, 'r','LineWidth',1); grid on; 
hold on
y2_PC = mean (sparsity_PC_Group2_50);
z2_PC = std (sparsity_PC_Group2_50)/sqrt (length (sparsity_PC_Group2_50));
errorbar (y2_PC,ci*z2_PC, 'b','LineWidth',1); grid on; 
legend('Pre TMS', 'Post TMS')
title('Absolute Participation Coefficient')
% % % figure;
% % % y1_GE = mean (Sparsity_GE_Group1_50);
% % % z1_GE = std (Sparsity_GE_Group1_50)/sqrt (length (Sparsity_GE_Group1_50)); 
% % % errorbar (y1_GE,ci*z1_GE, 'r','LineWidth',1); grid on; 
% % % hold on
% % % y2_GE = mean (Sparsity_GE_Group2_50);
% % % z2_GE = std (Sparsity_GE_Group2_50)/sqrt (length (Sparsity_GE_Group2_50));
% % % errorbar (y2_GE,ci*z2_GE, 'b','LineWidth',1); grid on; 
% % % legend('Pre TMS', 'Post TMS')
% % % title('Absolute Global Efficency')
% % % 
% % % figure;
% % % y1_BC = mean (sparsity_BC_Group1_50);
% % % z1_BC = std (sparsity_BC_Group1_50)/sqrt (length (sparsity_BC_Group1_50)); 
% % % errorbar (y1_BC,ci*z1_BC, 'r','LineWidth',1); grid on; 
% % % hold on
% % % y2_BC = mean (sparsity_BC_Group2_50);
% % % z2_BC = std (sparsity_BC_Group2_50)/sqrt (length (sparsity_BC_Group2_50));
% % % errorbar (y2_BC,ci*z2_BC, 'b','LineWidth',1); grid on; 
% % % legend('Pre TMS', 'Post TMS')
% % % title('Absolute Betweenness Centrality')
% % % 
% % % figure;
% % % y1_EC = mean (sparsity_EC_Group1_50);
% % % z1_EC = std (sparsity_EC_Group1_50)/sqrt (length (sparsity_EC_Group1_50)); 
% % % errorbar (y1_EC,ci*z1_EC, 'r','LineWidth',1); grid on; 
% % % hold on
% % % y2_EC = mean (sparsity_EC_Group2_50);
% % % z2_EC = std (sparsity_EC_Group2_50)/sqrt (length (sparsity_EC_Group2_50));
% % % errorbar (y2_EC,ci*z2_EC, 'b','LineWidth',1); grid on; 
% % % legend('Pre TMS', 'Post TMS')
% % % title('Absolute Eigenvector Centrality')
% % % 
% % % figure;
% % % y1_FC = mean (sparsity_FC_Group1_50);
% % % z1_FC = std (sparsity_FC_Group1_50)/sqrt (length (sparsity_FC_Group1_50)); 
% % % errorbar (y1_FC,ci*z1_FC, 'r','LineWidth',1); grid on; 
% % % hold on
% % % y2_FC = mean (sparsity_FC_Group2_50);
% % % z2_FC = std (sparsity_FC_Group2_50)/sqrt (length (sparsity_FC_Group2_50));
% % % errorbar (y2_FC,ci*z2_FC, 'b','LineWidth',1); grid on; 
% % % legend('Pre TMS', 'Post TMS')
% % % title('Absolute Flow Coefficient')
% % % 
% % % figure;
% % % y1_KC = mean (sparsity_KC_Group1_50);
% % % z1_KC = std (sparsity_KC_Group1_50)/sqrt (length (sparsity_KC_Group1_50)); 
% % % errorbar (y1_KC,ci*z1_KC, 'r','LineWidth',1); grid on; 
% % % hold on
% % % y2_KC = mean (sparsity_KC_Group2_50);
% % % z2_KC = std (sparsity_KC_Group2_50)/sqrt (length (sparsity_KC_Group2_50));
% % % errorbar (y2_KC,ci*z2_KC, 'b','LineWidth',1); grid on; 
% % % legend('Pre TMS', 'Post TMS')
% % % title('Absolute K-Core')

%% %% Random CC and PC ploting
% % figure;
% % y1_CC_rand = mean (sparsity_CC_rand_Group1_50);
% % z1_CC_rand = std (sparsity_CC_rand_Group1_50)/sqrt (length (sparsity_CC_rand_Group1_50)); 
% % errorbar (y1_CC_rand,ci*z1_CC_rand, 'r','LineWidth',1); grid on; 
% % hold on
% % y2_CC_rand = mean (sparsity_CC_rand_Group2_50);
% % z2_CC_rand = std (sparsity_CC_rand_Group2_50)/sqrt (length (sparsity_CC_rand_Group2_50));
% % errorbar (y2_CC_rand,ci*z2_CC_rand, 'b','LineWidth',1); grid on; 
% % title('Random Clustering Coefficient')
% % 
% % figure;
% % y1_PC_rand = mean (sparsity_PC_rand_Group1_50);
% % z1_PC_rand = std (sparsity_PC_rand_Group1_50)/sqrt (length (sparsity_PC_rand_Group1_50)); 
% % errorbar (y1_PC_rand,ci*z1_PC_rand, 'r','LineWidth',1); grid on; 
% % hold on
% % y2_PC_rand = mean (sparsity_PC_rand_Group2_50);
% % z2_PC_rand = std (sparsity_PC_rand_Group2_50)/sqrt (length (sparsity_PC_rand_Group2_50));
% % errorbar (y2_PC_rand,ci*z2_PC_rand, 'b','LineWidth',1); grid on; 
% % title('Random Participation Coefficient')
% % % figure; plot(sparsity_PC_Group1_50')
% % % figure; plot(sparsity_PC_Group2_50')
% % % figure; plot(sparsity_PC_rand_Group1_50')
% % % figure; plot(sparsity_PC_rand_Group2_50')
%% %Ploting Normalized CC, PC, PL, GE, LE Images

figure;
y1_CC = mean (sparsity_CC_normalised_Group1);
z1_CC = std (sparsity_CC_normalised_Group1)/sqrt (length (sparsity_CC_normalised_Group1)); 
errorbar (y1_CC,ci*z1_CC, 'r','LineWidth',1); grid on; 
hold on
y2_CC = mean (sparsity_CC_normalised_Group2);
z2_CC = std (sparsity_CC_normalised_Group2)/sqrt (length (sparsity_CC_normalised_Group2));
errorbar (y2_CC,ci*z2_CC, 'b','LineWidth',1); grid on; 
legend('Pre TMS', 'Post TMS')
title('Normalised Clustering Coefficient')
xticks([1 3 5 7 9 11 13 15]); xticklabels({'0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.55'})
set(gca,'FontSize',12,'FontWeight','B')
%%%
figure;
hold on
y1_PL = mean (sparsity_PL_normalised_Group1);
z1_PL = std (sparsity_PL_normalised_Group1)/sqrt (length (sparsity_PL_normalised_Group1)); 
errorbar (y1_PL,ci*z1_PL, 'r','LineWidth',1); grid on; 
hold on
y2_PL = mean (sparsity_PL_normalised_Group2);
z2_PL = std (sparsity_PL_normalised_Group2)/sqrt (length (sparsity_PL_normalised_Group2));
errorbar (y2_PL,ci*z2_PL, 'b','LineWidth',1); grid on;
legend('Pre TMS', 'Post TMS')
title('Normalised Path Length')
xticks([1 3 5 7 9 11 13 15]); xticklabels({'0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.55'})
set(gca,'FontSize',12,'FontWeight','B')

figure;
y1_SW = mean (SmallWorldNess_Group1);
z1_SW = std (SmallWorldNess_Group1)/sqrt (length (SmallWorldNess_Group1));
errorbar (y1_SW,ci*z1_SW, 'r','LineWidth',1); grid on;
hold on
y2_SW = mean (SmallWorldNess_Group2);
z2_SW = std (SmallWorldNess_Group2)/sqrt (length (SmallWorldNess_Group2));
errorbar (y2_SW,ci*z2_SW, 'b','LineWidth',1); grid on;
legend('Pre TMS', 'Post TMS')
title('Small Worldness')
xticks([1 3 5 7 9 11 13 15]); xticklabels({'0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.55'})
set(gca,'FontSize',12,'FontWeight','B')

figure;
y1_CC = mean (sparsity_LE_normalised_Group1);
z1_CC = std (sparsity_LE_normalised_Group1)/sqrt (length (sparsity_LE_normalised_Group1)); 
errorbar (y1_CC,ci*z1_CC, 'r','LineWidth',1); grid on; 
hold on
y2_CC = mean (sparsity_LE_normalised_Group2);
z2_CC = std (sparsity_LE_normalised_Group2)/sqrt (length (sparsity_LE_normalised_Group2));
errorbar (y2_CC,ci*z2_CC, 'b','LineWidth',1); grid on; 
legend('Pre TMS', 'Post TMS')
title('Normalised Local Eficency')
xticks([1 3 5 7 9 11 13 15]); xticklabels({'0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.55'})
set(gca,'FontSize',12,'FontWeight','B')

figure;
hold on
y1_PL = mean (sparsity_GE_normalised_Group1);
z1_PL = std (sparsity_GE_normalised_Group1)/sqrt (length (sparsity_GE_normalised_Group1)); 
errorbar (y1_PL,ci*z1_PL, 'r','LineWidth',1); grid on; 
hold on
y2_PL = mean (sparsity_GE_normalised_Group2);
z2_PL = std (sparsity_GE_normalised_Group2)/sqrt (length (sparsity_GE_normalised_Group2));
errorbar (y2_PL,ci*z2_PL, 'b','LineWidth',1); grid on;
legend('Pre TMS', 'Post TMS')
title('Normalised Global Efficency')
xticks([1 3 5 7 9 11 13 15]); xticklabels({'0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.55'})
set(gca,'FontSize',12,'FontWeight','B')
%%%
figure;
y1_PC = mean (sparsity_PC_normalised_Group1);
z1_PC = std (sparsity_PC_normalised_Group1)/sqrt (length (sparsity_PC_normalised_Group1)); 
errorbar (y1_PC,ci*z1_PC, 'r','LineWidth',1); grid on;  
hold on
y2_PC = mean (sparsity_PC_normalised_Group2);
z2_PC = std (sparsity_PC_normalised_Group2)/sqrt (length (sparsity_PC_normalised_Group2));
errorbar (y2_PC,ci*z2_PC, 'b','LineWidth',1); grid on; 
legend('Pre TMS', 'Post TMS')
title('Normalised Participation Coefficient')
xticks([1 3 5 7 9 11 13 15]); xticklabels({'0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.55'})
set(gca,'FontSize',12,'FontWeight','B')

%%
figure;
sub1=1:4; sub2=5:8; sub3=9:12; sub4=13:16; sub5=17:20; sub6=21:24;
sub=sub2;
plot(sparsity_GE_normalised_Group1(sub,:)','b'); 
hold on
plot(sparsity_GE_normalised_Group2(sub,:)','r');
%% % Plot Normalised Betweenness Centrality
% % % figure;
% % % y1_BC = mean (sparsity_BC_normalised_Group1);
% % % z1_BC = std (sparsity_BC_normalised_Group1)/sqrt (length (sparsity_BC_normalised_Group1)); 
% % % errorbar (y1_BC,ci*z1_BC, 'r','LineWidth',1); grid on;  
% % % hold on
% % % y2_BC = mean (sparsity_BC_normalised_Group2);
% % % z2_BC = std (sparsity_BC_normalised_Group2)/sqrt (length (sparsity_BC_normalised_Group2));
% % % errorbar (y2_BC,ci*z2_BC, 'b','LineWidth',1); grid on; 
% % % legend('Pre TMS', 'Post TMS')
% % % title('Normalised Betweenness Centrality')
% % % xticks([1 3 5 7 9 11 13 15]); xticklabels({'0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.55'})
% % % set(gca,'FontSize',12,'FontWeight','B')

% % % %%% Plot Normalised Eigenvector Centrality'
% % % figure;
% % % y1_EC = mean (sparsity_EC_normalised_Group1);
% % % z1_EC = std (sparsity_EC_normalised_Group1)/sqrt (length (sparsity_EC_normalised_Group1)); 
% % % errorbar (y1_EC,ci*z1_EC, 'r','LineWidth',1); grid on;  
% % % hold on
% % % y2_EC = mean (sparsity_EC_normalised_Group2);
% % % z2_EC = std (sparsity_EC_normalised_Group2)/sqrt (length (sparsity_EC_normalised_Group2));
% % % errorbar (y2_EC,ci*z2_EC, 'b','LineWidth',1); grid on; 
% % % legend('Pre TMS', 'Post TMS')
% % % title('Normalised Eigenvector Centrality')
% % % xticks([1 3 5 7 9 11 13 15]); xticklabels({'0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.55'})
% % % set(gca,'FontSize',12,'FontWeight','B')

%%% Plot Normalised Flow Coefficient
% % % figure;
% % % y1_FC = mean (sparsity_FC_normalised_Group1);
% % % z1_FC = std (sparsity_FC_normalised_Group1)/sqrt (length (sparsity_FC_normalised_Group1)); 
% % % errorbar (y1_FC,ci*z1_FC, 'r','LineWidth',1); grid on;  
% % % hold on
% % % y2_FC = mean (sparsity_FC_normalised_Group2);
% % % z2_FC = std (sparsity_FC_normalised_Group2)/sqrt (length (sparsity_FC_normalised_Group2));
% % % errorbar (y2_FC,ci*z2_FC, 'b','LineWidth',1); grid on; 
% % % legend('Pre TMS', 'Post TMS')
% % % title('Normalised Flow Coefficient')
% % % xticks([1 3 5 7 9 11 13 15]); xticklabels({'0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.55'})
% % % set(gca,'FontSize',12,'FontWeight','B')

% % % Plot Normalised K Core
% % % figure;
% % % y1_KC = mean (sparsity_KC_normalised_Group1);
% % % z1_KC = std (sparsity_KC_normalised_Group1)/sqrt (length (sparsity_KC_normalised_Group1)); 
% % % errorbar (y1_KC,ci*z1_KC, 'r','LineWidth',1); grid on;  
% % % hold on
% % % y2_KC = mean (sparsity_KC_normalised_Group2);
% % % z2_KC = std (sparsity_KC_normalised_Group2)/sqrt (length (sparsity_KC_normalised_Group2));
% % % errorbar (y2_KC,ci*z2_KC, 'b','LineWidth',1); grid on; 
% % % legend('Pre TMS', 'Post TMS')
% % % title('Normalised K Core')
% % % xticks([1 3 5 7 9 11 13 15]); xticklabels({'0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.55'})
% % % set(gca,'FontSize',12,'FontWeight','B')
%%% Plot Modularity
% % % figure;
% % % hold on
% % % y1_Mod = mean (sparsity_Mod_normalised_Group1);
% % % z1_Mod = std (sparsity_Mod_normalised_Group1)/sqrt (length (sparsity_Mod_normalised_Group1)); 
% % % errorbar (y1_Mod,ci*z1_Mod, 'r','LineWidth',1); grid on; 
% % % hold on
% % % y2_Mod = mean (sparsity_Mod_normalised_Group2);
% % % z2_Mod = std (sparsity_Mod_normalised_Group2)/sqrt (length (sparsity_Mod_normalised_Group2));
% % % errorbar (y2_Mod,ci*z2_Mod, 'b','LineWidth',1); grid on;
% % % legend('Pre TMS', 'Post TMS')
% % % title('Normalised Modularity')
% % % xticks([1 3 5 7 9 11 13 15]); xticklabels({'0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.55'})
% % % set(gca,'FontSize',12,'FontWeight','B')
% % % %%  Plot Normalised Assortativity
% % % figure;
% % % hold on
% % % y1_AS = mean (sparsity_AS_normalised_Group1);
% % % z1_AS = std (sparsity_AS_normalised_Group1)/sqrt (length (sparsity_AS_normalised_Group1)); 
% % % errorbar (y1_AS,ci*z1_AS, 'r','LineWidth',1,'LineWidth',1); grid on; 
% % % hold on
% % % y2_AS = mean (sparsity_AS_normalised_Group2);
% % % z2_AS = std (sparsity_AS_normalised_Group2)/sqrt (length (sparsity_AS_normalised_Group2));
% % % errorbar (y2_AS,ci*z2_AS, 'b','LineWidth',1,'LineWidth',1); grid on;
% % % legend('Pre TMS', 'Post TMS')
% % % title('Normalised Assortativity')
% % % xticks([1 3 5 7 9 11 13 15]); xticklabels({'0.9','0.85','0.8','0.75','0.7','0.65','0.6','0.55'})
% % % set(gca,'FontSize',12,'FontWeight','B')

%%
Grp_con = ones(length(SUBJlist_Group1),1); Grp_mcs = 2*(ones(length(SUBJlist_Group2),1)); 
Grp = [Grp_con; Grp_mcs]; 
Grp_PC = [mean(sparsity_PC_normalised_Group1,2); mean(sparsity_PC_normalised_Group2,2)];
figure; notBoxPlot(Grp_PC,Grp,0.5,'patch',ones(length(Grp_PC),1)); title('PC')
Grp_CC = [mean(sparsity_CC_normalised_Group1,2); mean(sparsity_CC_normalised_Group2,2)];
figure; notBoxPlot(Grp_CC,Grp,0.5,'patch',ones(length(Grp_CC),1)); title('CC')
Grp_GE = [mean(sparsity_GE_normalised_Group1,2); mean(sparsity_GE_normalised_Group2,2)]; 
figure; notBoxPlot(Grp_GE,Grp,0.5,'patch',ones(length(Grp_GE),1)); title('GE')
Grp_LE = [mean(sparsity_LE_normalised_Group1,2); mean(sparsity_LE_normalised_Group2,2)]; 
figure; notBoxPlot(Grp_LE,Grp,0.5,'patch',ones(length(Grp_LE),1)); title('LE')
Grp_con = ones(length(SUBJlist_Group1),1); Grp_mcs = 2*(ones(length(SUBJlist_Group2),1)); 
Grp_EC = [mean(sparsity_EC_normalised_Group1,2); mean(sparsity_EC_normalised_Group2,2)];
figure; notBoxPlot(Grp_EC,Grp,0.5,'patch',ones(length(Grp_EC),1)); title('EC')
%%
for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_NISI_normalised_Group1(i,j) = (sparsity_PC_normalised_Group1(i,j).*sparsity_GE_normalised_Group1(i,j))/(sparsity_CC_normalised_Group1(i,j).*sparsity_LE_normalised_Group1(i,j)); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_NISI_normalised_Group2(i,j) = (sparsity_PC_normalised_Group2(i,j).*sparsity_GE_normalised_Group2(i,j))/(sparsity_CC_normalised_Group2(i,j).*sparsity_LE_normalised_Group2(i,j)); 
    end
end
%%
for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_NISI_normalised_Group1(i,j) = (sparsity_PC_normalised_Group1(i,j))/(sparsity_CC_normalised_Group1(i,j)); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_NISI_normalised_Group2(i,j) = (sparsity_PC_normalised_Group2(i,j))/(sparsity_CC_normalised_Group2(i,j)); 
    end
end
figure;
y1_PC = mean (sparsity_NISI_normalised_Group1);
z1_PC = std (sparsity_NISI_normalised_Group1)/sqrt (length (sparsity_NISI_normalised_Group1)); 
errorbar (y1_PC,ci*z1_PC, 'r','LineWidth',1); grid on;  
hold on
y2_PC = mean (sparsity_NISI_normalised_Group2);
z2_PC = std (sparsity_NISI_normalised_Group2)/sqrt (length (sparsity_NISI_normalised_Group2));
errorbar (y2_PC,ci*z2_PC, 'b','LineWidth',1); grid on; 
title('Integration Segregation Index')

Grp_NISI = [squeeze(mean(sparsity_NISI_normalised_Group1,2)); squeeze(mean(sparsity_NISI_normalised_Group2,2))];
figure; notBoxPlot(Grp_NISI,Grp,0.5,'patch',ones(length(Grp_NISI),1));
%%
for i = 1:length(SUBJlist_Group1); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_NISI1_normalised_Group1(i,j) = sparsity_PC_normalised_Group1(i,j)/sparsity_LE_normalised_Group1(i,j); 
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = SparcityStart:SparcityEnd; 
            sparsity_NISI1_normalised_Group2(i,j) = sparsity_PC_normalised_Group2(i,j)/sparsity_LE_normalised_Group2(i,j); 
    end
end

figure;
y1_PC = mean (sparsity_NISI1_normalised_Group1);
z1_PC = std (sparsity_NISI1_normalised_Group1)/sqrt (length (sparsity_NISI1_normalised_Group1)); 
errorbar (y1_PC,ci*z1_PC, 'r','LineWidth',1); grid on;  
hold on
y2_PC = mean (sparsity_NISI1_normalised_Group2);
z2_PC = std (sparsity_NISI1_normalised_Group2)/sqrt (length (sparsity_NISI1_normalised_Group2));
errorbar (y2_PC,ci*z2_PC, 'b','LineWidth',1); grid on; 
title('Integration Segregation Index')

Grp_NISI = [squeeze(mean(sparsity_NISI1_normalised_Group1,2)); squeeze(mean(sparsity_NISI1_normalised_Group2,2))];
figure; notBoxPlot(Grp_NISI,Grp,0.5,'patch',ones(length(Grp_NISI),1));
%%
for i = SparcityStart:SparcityEnd
    [h_CC1(i),p_CC1(i)] = ttest2(sparsity_CC_normalised_Group2(:,i),sparsity_CC_normalised_Group1(:,i),0.05,'left');
end

h_CC1
p_CC1
%%
for i = SparcityStart:SparcityEnd
    [h_SW1(i),p_SW1(i)] = ttest2(SmallWorldNess_Group2(:,i),SmallWorldNess_Group1(:,i),0.05,'left');
end

h_SW1
p_SW1
%%
for i = SparcityStart:SparcityEnd
    [h_PC1(i),p_PC1(i)] = ttest2(sparsity_PC_normalised_Group2(:,i),sparsity_PC_normalised_Group1(:,i),0.05,'right');
end

h_PC1
p_PC1
%%
%----------------------%% Brain resion significant Computations for CC %-------------%
clc; spr=logical([0 1 1 1 1 1 1 1 0]);

sparsity_CC_Group1_ROI = squeeze(mean(Group1_CC_50(:,spr,:),2)); 
sparsity_CC_Group2_ROI = squeeze(mean(Group2_CC_50(:,spr,:),2));
sparsity_CC_rand_Group1_ROI = squeeze(mean(Group1_CC_rand_squ(:,spr,:),2));
sparsity_CC_rand_Group2_ROI = squeeze(mean(Group2_CC_rand_squ(:,spr,:),2));

for i = 1:length(SUBJlist_Group1); 
    for j = 1:NoOfRoi; 
            sparsity_CC_normalised_Group1_ROI(i,j) = sparsity_CC_Group1_ROI(i,j)/sparsity_CC_rand_Group1_ROI(i,j);
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:NoOfRoi; 
            sparsity_CC_normalised_Group2_ROI(i,j) = sparsity_CC_Group2_ROI(i,j)/sparsity_CC_rand_Group2_ROI(i,j); 
    end
end

for i = 1:NoOfRoi
    [h_CC_ROI_Normalised(i),p_CC_ROI_Normalised(i)] = ttest2(sparsity_CC_normalised_Group1_ROI(:,i),sparsity_CC_normalised_Group2_ROI(:,i),0.01,'right');
end

h_CC_ROI_Normalised
p_CC_ROI_Normalised
%%
%----------------------%% Brain region significant Computations for PC %-------------%
clc; spr=logical([0 1 1 1 1 1 1 1 0]);

sparsity_PC_Group1_ROI = squeeze(mean(Group1_PC_50(:,spr,:),2)); 
sparsity_PC_Group2_ROI = squeeze(mean(Group2_PC_50(:,spr,:),2)); 
sparsity_PC_rand_Group1_ROI = squeeze(mean(Group1_PC_rand_squ(:,spr,:),2));
sparsity_PC_rand_Group2_ROI = squeeze(mean(Group2_PC_rand_squ(:,spr,:),2));

for i = 1:length(SUBJlist_Group1); 
    for j = 1:NoOfRoi; 
            sparsity_PC_normalised_Group1_ROI(i,j) = sparsity_PC_Group1_ROI(i,j)/sparsity_PC_rand_Group1_ROI(i,j);
    end
end

for i = 1:length(SUBJlist_Group2); 
    for j = 1:NoOfRoi; 
            sparsity_PC_normalised_Group2_ROI(i,j) = sparsity_PC_Group2_ROI(i,j)/sparsity_PC_rand_Group2_ROI(i,j); 
    end
end

for i = 1:NoOfRoi
    [h_PC_ROI_Normalised(i),p_PC_ROI_Normalised(i)] = ttest2(sparsity_PC_normalised_Group2_ROI(:,i),sparsity_PC_normalised_Group1_ROI(:,i),0.01,'right');
end

h_PC_ROI_Normalised   
p_PC_ROI_Normalised
%%[p,h]=fdr(p_CC_ROI_Normalised,0.05);
%p

%%  Graphical Plot of GT measures %%%
% % % tpz_cc_control = trapz(sparsity_CC_normalised_Control(:,1:30),2)/30;
% % % tpz_cc_Patient = trapz(sparsity_CC_normalised_Patient(:,1:30),2)/30;
% % % tpz_pl_control = trapz(sparsity_PL_normalised_Control(:,1:30),2)/30;
% % % tpz_pl_Patient = trapz(sparsity_PL_normalised_Patient(:,1:30),2)/30;
% % % tpz_sw_control = trapz(SmallWorldNess_Control(:,1:30),2)/30;
% % % tpz_sw_Patient = trapz(SmallWorldNess_Patient(:,1:30),2)/30;
% % % tpz_le_control = trapz(sparsity_LE_normalised_Control(:,1:30),2)/30;
% % % tpz_le_Patient = trapz(sparsity_LE_normalised_Patient(:,1:30),2)/30;
% % % tpz_ge_control = trapz(sparsity_GE_normalised_Control(:,1:30),2)/30;
% % % tpz_ge_Patient = trapz(sparsity_GE_normalised_Patient(:,1:30),2)/30;
% % % 
% % % tpz_GT = [tpz_cc_control tpz_cc_Patient tpz_pl_control tpz_pl_Patient tpz_sw_control tpz_sw_Patient tpz_le_control tpz_le_Patient tpz_ge_control tpz_ge_Patient]
% % % 
% % % [p,h]=ttest2(tpz_GT(:,1),tpz_GT(:,2), 0.05,'left')