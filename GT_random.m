%%%%%%%%%%%%%% To compute graph theory measures of FC matrix %%%%%%%%%%%%
%%%%%%% Author: Sneha Ray, N3 Lab, UCSF, US
%%%%%%% Supervised: Prof. Andrew Moses Lee
%%%%%%% Last Date of Modification: 1 Dec 2024
%%%Software/Toolbox- BCT and Matlab
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc ]
%%% Specify path of data folder and group name
path='/Users/snray/Documents/data_lab/AN/BoldData/';
cd(path)
SUBJlist=dir('AN*');
group ='post';
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
for i=1                                                                                            :length(SUBJlist)
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
        sparsity_val=0.125:0.025:0.475;% 0.1:0.05:0.9; %sparsity_val=0.01:0.025:0.52; %sparsity_val=0.01:0.025:1;
        for i2=1:length (sparsity_val)                   %i2=1:40
            %% %%%Network properties/ network measurement 
            GT_sparsity(i2)=sparsity_val(i2); % define sparsity
            corr_data_thr1=threshold_proportional(GT_corr_data_abs,GT_sparsity(i2)); % calcualte binary matrix
            corr_data_thr_bin=weight_conversion(corr_data_thr1,'binarize'); % binary weight matrix
            corr_data_thr = weight_conversion(corr_data_thr_bin, 'autofix'); % removing NaN & Inf
                       
            for random_number=1:10   % random_number=1:50
                             
                random_network=randmio_und(corr_data_thr,5);
                %%
                GT_corr_data_rand_thr(i2,random_number,:,:)=random_network;% asign to different array
                
                GT_degree_rand(i2,random_number,:)=degrees_und(random_network); % calculate degree
                
                GT_clust_coeff_rand(i2,random_number,:)=clustering_coef_bu(random_network);%% calculate clustering coeff
                                
                GT_local_eff_rand(i2,random_number,:)=efficiency_bin(random_network,1); % local efficiency
                
                GT_global_eff_rand(i2,random_number,:)=efficiency_bin(random_network,0); % global efficiency
                
                GT_distance_matrix_rand(i2,random_number,:,:)=distance_bin(random_network); % distance matrix
                
                GT_path_length_rand(i2,random_number) =charpath(squeeze(GT_distance_matrix_rand(i2,random_number,:,:)),1,0); % path length

                GT_betw_cent_rand(i2,random_number,:)= betweenness_bin(random_network); % betweenness centrality
                        
                GT_eigv_cent_rand(i2,random_number,:)= eigenvector_centrality_und(random_network); % Eigenvector centrality
            
                GT_assortvity_rand (i2,random_number) = assortativity_bin(random_network,0); % assortativity
            
                GT_flow_coef_rand(i2,random_number,:)= flow_coef_bd(random_network); % flow coefficient is similar to betweenness centrality

                GT_k_core_rand(i2,random_number,:)= kcoreness_centrality_bu(random_network); % k_core
            
            
            %%%%% Participation coefficient and modspan
                param.heuristic=50;
                for i = 1:param.heuristic
                    [Ci, allQ(i2,random_number,i)] = community_louvain(random_network);
                     allCi(i2,random_number,i,:) = Ci;
                     allpc(i2,random_number,i,:) = participation_coef(random_network,Ci); 
                end
                GT_modularity_rand(i2,random_number)= mean(allQ(i2,random_number,:));  % modularity
                GT_community_structure_rand(i2,random_number,1:ROIs) = squeeze(allCi(i2,random_number,1,:)); % community structure
                GT_participation_coeff_rand(i2,random_number,1:ROIs) = mean(squeeze(allpc(i2,random_number,:,:)));  %participation coefficient
                %%%%%    
            end
        end
%%
        SUBJname1=  SUBJname;
        varname=([SUBJname '_' group '_RAND'])
        save(varname);
        cd ..
        cd ..
end
