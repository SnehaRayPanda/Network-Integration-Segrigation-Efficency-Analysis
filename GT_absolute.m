%%%%%%%%%%%%%% To compute graph theory measures of FC matrix %%%%%%%%%%%%
%%%%%%% Author: Sneha Ray, N3 Lab, UCSF, US
%%%%%%% Supervised: Prof. Andrew Moses Lee
%%%%%%% Last Date of Modification: 1 Dec 2024
%%%Software/Toolbox- BCT and Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
%%% Specify path of data folder and group name
path='/Users/snray/Documents/data_lab/AN/BoldData/';
cd(path)
SUBJlist=dir('AN*');
group ='pre';
%SUBJlist=SUBJlist(1:length(SUBJlist));
%% filter Initialisation
     flp = 0.03;   %flp = .01;  0.03      % lowpass frequency of filter
     fhi = 0.07;   %fhi = .09;  0.07       % highpass
     delta=0.8;    %TR
     k = 2;                  % 2nd order butterworth filter
     fnq = 1/(2*delta);       % Nyquist frequency
     Wn = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
     [bfilt2,afilt2] = butter(k,Wn); 
%%
for i=1:length(SUBJlist)
    SUBJname=SUBJlist(i).name
    cd(SUBJname)
    cd(group)
        %%% load the data
        File_ID=dir('*tsv');
        File_ID=File_ID.name;

        ts_4s456_t = readtable(File_ID, 'FileType', 'text', 'Delimiter', '\t');
        GT_data= table2array(ts_4s456_t); % BOLD time series

         for seed = 1:size(GT_data,2)
          x = GT_data(:,seed);
          GT_data(:,seed) = filtfilt(bfilt2,afilt2,x);
         end 
         
        GT_corr_data=corr(GT_data); % create NxN connectivity matrix (N= No of ROIs/Channels)
        GT_corr_data_abs = abs(GT_corr_data); % conver all connectivity values to absolute values
        ROIs = size(GT_corr_data_abs,2); % No  of ROIs      
        %% Thresholding
        sparsity_val=0.1:0.05:0.5;% 0.1:0.05:0.9; %sparsity_val=0.01:0.025:0.52; %sparsity_val=0.01:0.025:1;
        for i2=1:length (sparsity_val)                   %i2=1:40
            %% %%%Network properties / network measurement 
            GT_sparsity(i2)=sparsity_val(i2); % define sparsity
            corr_data_thr1=threshold_proportional(GT_corr_data_abs,GT_sparsity(i2)); % compute the thresholded matrix
            corr_data_thr_bin=weight_conversion(corr_data_thr1,'binarize'); % convert binary weight matrix
            corr_data_thr = weight_conversion(corr_data_thr_bin, 'autofix'); % removing NaN & Inf
            %%
            GT_corr_data_thr(i2,:,:)=corr_data_thr;% asign to different array
            
            GT_degree(i2,:)=degrees_und(corr_data_thr); % calculate degree
            
            GT_clust_coeff(i2,:)=clustering_coef_bu(corr_data_thr);%% calculate clustering coeff
            
            GT_local_eff(i2,:)=efficiency_bin(corr_data_thr,1); % local efficiency
            
            GT_global_eff(i2,:)=efficiency_bin(corr_data_thr,0); % global efficiency
            
            GT_distance_matrix(i2,:,:)=distance_bin(corr_data_thr); % distance matrix
                       
            GT_path_length(i2)=charpath(squeeze(GT_distance_matrix(i2,:,:)),1,0); % path length
            
            GT_betw_cent(i2,:)= betweenness_bin(corr_data_thr); % betweenness centrality
                        
            GT_eigv_cent(i2,:)= eigenvector_centrality_und(corr_data_thr); % Eigenvector centrality
            
            GT_assortvity (i2) = assortativity_bin(corr_data_thr,0); % assortativity
            
            GT_flow_coef(i2,:)= flow_coef_bd(corr_data_thr); % flow coefficient is similar to betweenness centrality

            GT_k_core(i2,:)= kcoreness_centrality_bu(corr_data_thr); % k_core
                      
              
            %%%%% Participation coefficient and modspan
            param.heuristic=50;
            
            for i = 1:param.heuristic
                [Ci, allQ(i2,i)] = community_louvain(corr_data_thr);
            
                allCi(i2,i,:) = Ci;
   
                allpc(i2,i,:) = participation_coef(corr_data_thr,Ci); 
            end
        
            GT_modularity(i2)= mean(allQ(i2,:));  % modularity
            GT_community_structure(i2,1:ROIs) = squeeze(allCi(i2,1,:)); % community structure
            GT_participation_coeff(i2,1:ROIs) = mean(squeeze(allpc(i2,:,:)));  %participation coefficient
            %%%%%
            
        end
 %%     
        SUBJname1=  SUBJname;
        varname=([SUBJname '_' group '_ABS'])
        save(varname);
        cd ..
        cd ..
end