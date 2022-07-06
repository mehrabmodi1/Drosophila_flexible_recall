clear all
close all

dataset_list_paths = [ ...
                         %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\d5HT1b_similar_od_handovers.xls'};...
                         %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\c305a_similar_od_handovers.xls'};...
                         {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\c739_similar_od_handovers.xls'};...
                      ];

                  
                  
[del, odor_names1] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
[del, odor_names2] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);
paper_save_dir = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\Fig2_KCs_simple\';
paper_save_dir_tr = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\Fig5_KCs_transitions\';
paper_save_dir_sfig = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\SFig1_KCs_simple_analyses\';
odor_names2{3} = 'Butyl acetate';
od_names = [{'PA'}, {'BA'}, {'EL'}];
           
paired_color = [0,136,55]./256;
unpaired_color = [166,219,160]./256;
EL_color = [123,50,148]./256;
PA_color = [0,136,55]./256;
BA_color = [166,219,160]./256;
mean_color = [0.8, 0.4, 0.4];


suppress_QC_plots = 1;
pool_flies_traces = 0;

global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;

y_ax_traces = 0.8;
y_ax_fit_traces = 0.6;
saved_long_traces = 0;
all_sig_frs = [];
pause_PCAs = 0;
single_fly_n = 5;   %G KCs - 5, A'\B' KCs - 5, A/B KCs - 3
subtract_pulse1_resps = 0;      %1 - subtracts the appropriate mean responses to single pulses from the transition response traces to look for linear summing

if isempty(findstr(dataset_list_paths{1}, 'd5HT1b')) == 0
    KC_type = 'Gamma';
elseif isempty(findstr(dataset_list_paths{1}, 'c739')) == 0
    KC_type = 'AlphaBeta';
elseif isempty(findstr(dataset_list_paths{1}, 'c305a')) == 0
    KC_type = 'AlphapBetap';
else
end

for list_n = 1:size(dataset_list_paths, 1)
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    curr_dir_list_path = update_list_path(curr_dir_list_path);
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'list_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));
        
    %loop to go through all experiment datasets listed in list file
    all_resp_traces = [];
    sig_cell_mat_all = [];
    saved_resp_sizes_all = [];
    n_cells_all = [];
    all_fly_accuracy_mat = [];
    for dir_n = 1:n_dirs
        fly_n = fly_n + 1;
              
        saved_an_results.scriptname = mfilename('fullpath');
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, 3);
        curr_dir = [curr_dir, '\1\'];
       
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        
        cd(curr_dir);
        tif_name = dir('*.tif');
        stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
        [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);
        [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig, PID_traces_matched, PID_traces_orig] = load_params_trains_modular(curr_dir, tif_times, frame_time);    %reading in trial stimulus parameters after matching time stamps to F traces
        color_vec = [PA_color; BA_color; EL_color];
        
        %Reading in experimental parameters
        odor_list_olf1 = unique(stim_mat_simple(:, 1) );
        n_odors_olf1 = length(odor_list_olf1);
        odor_dur_list_olf1 = unique(stim_mat_simple(:, 2) );
        n_od_durs_olf1 = length(odor_dur_list_olf1);
        n_trains_olf1 = max(stim_mat_simple(:, 8));
        
        odor_list_olf2 = unique(stim_mat_simple(:, 9) );
        n_odors_olf2 = length(odor_list_olf2);
        odor_dur_list_olf2 = unique(stim_mat_simple(:, 10) );
        n_od_durs_olf2 = length(odor_dur_list_olf2);
        n_trains_olf2 = max(stim_mat_simple(:, 17));
        
        
        
        %loading extracted raw fluorescence data matrices written by raw_dff_extractor
        raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        tif_n_col_n = find_stim_mat_simple_col('matched_tif_n', column_heads);
        
        %inserting dummy trials (all nans) into raw_data_mat for pairing trials for
        %which no corress .tifs were acquired
        raw_data_mat = match_up_rawmat_matchedtrs(raw_data_mat, stim_mat_simple, tif_n_col_n);
        
        %dumping data from manually identified, z-drifted trials
        bad_tr_list = 1:size(raw_data_mat, 3);
        bad_tr_list(good_tr_list) = [];
        raw_data_mat(:, :, bad_tr_list) = nan;
        
        n_cells = size(raw_data_mat, 2);
        
        %calculating dF/F traces from raw data
        filt_time = .5;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir, []);
        
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
        
        
        %identifying stim_mat_simple col numbers
        led_on_col_n = find_stim_mat_simple_col('led_on', column_heads);            %identifying relevant column number in stim_mat_simple
        od_olf1_col_n = find_stim_mat_simple_col('odor_n', column_heads);           %identifying relevant column number in stim_mat_simple
        od_olf2_col_n = find_stim_mat_simple_col('odour_olf2', column_heads);       %identifying relevant column number in stim_mat_simple
        dur_olf1_col_n = find_stim_mat_simple_col('duration', column_heads);        %identifying relevant column number in stim_mat_simple
        dur_olf2_col_n = find_stim_mat_simple_col('duration_olf2', column_heads);   %identifying relevant column number in stim_mat_simple
        od_col_ns = [od_olf1_col_n, od_olf2_col_n];
        dur_col_ns = [dur_olf1_col_n, dur_olf2_col_n];
        od_durs = unique(stim_mat_simple(:, dur_col_ns(2)));
        od_durs(isnan(od_durs)) = [];
        odn_list_olf2 = unique(stim_mat_simple(:, od_col_ns(2)));
        odn_list_olf2(isnan(odn_list_olf2)) = [];
        odn_list_olf1 = unique(stim_mat_simple(:, od_col_ns(1)));
        
        %replacing Nans in stim_mat_simple with zeroes
        stim_mat_simple_orig = stim_mat_simple;
        del = isnan(stim_mat_simple);
        stim_mat_simple(del) = 0;
        
        %replacing small vals in olf1_dur column with zeroes
        del = find(stim_mat_simple(:, dur_col_ns(1)) < 1);
        stim_mat_simple(del, dur_col_ns(1)) = 0;
       
        [resp_sizes, sig_trace_mat, sig_cell_mat, sig_cell_mat_key, resp_areaundercurves] = cal_sig_responses_res_modular(dff_data_mat, stim_mat, stim_mat_simple, frame_time, od_col_ns, dur_col_ns);
        
        %sanity-checking KC data
        [del, all_bad_trs] = cell_data_quality_control(dff_data_mat_f, stim_mat, stim_mat_simple, column_heads, sig_cell_mat, sig_cell_mat_key, suppress_QC_plots, frame_time);
                
        union_sig_cell_vec = sum(sig_cell_mat, 2);
        union_sig_cell_vec = sign(union_sig_cell_vec);
        union_sig_cells = find(union_sig_cell_vec == 1);
        sum(sig_cell_mat)./size(sig_cell_mat, 1);
        
        sig_cell_mat_all = [sig_cell_mat_all; sig_cell_mat(:, 1:3)];
        
        %computing mean response trace for each simple odor stimulus presented
        resp_trace_mat = [];
        bad_cells_all = [];
        saved_resp_sizes = [];
        single_traces_mat = [];
        n_reps_vec = [];
        for stim_type_n = 1:size(sig_cell_mat_key, 1)
            curr_stim_vec = sig_cell_mat_key(stim_type_n, :);
            
            olf1_od_n = curr_stim_vec(1, 1);
            olf1_dur = curr_stim_vec(1, 2);
            olf2_od_n = curr_stim_vec(1, 3);
            olf2_dur = curr_stim_vec(1, 4);
            
            if olf1_dur ~= 0 && olf2_dur ~= 0
                curr_stim_type = 2;     %handover stimuli
            elseif olf2_dur == 0
                curr_stim_type = 0;     %simple stimulus, delivered on olf1
                curr_stim_type_simp = 0;
            elseif olf1_dur == 0
                curr_stim_type = 1;     %simple stimulus, delivered on olf2
                curr_stim_type_simp = 1;
            end
            
            curr_sig_cells = find(sig_cell_mat(:, stim_type_n) == 1);
            curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == olf1_od_n & stim_mat_simple(:, dur_col_ns(1)) == olf1_dur &... 
                                stim_mat_simple(:, od_col_ns(2)) == olf2_od_n & stim_mat_simple(:, dur_col_ns(2)) == olf2_dur);
            
            %logging traces sorted by stim type to use for logistic regression decoding.
            single_traces_mat = pad_n_concatenate(single_traces_mat, dff_data_mat_f(:, :, curr_trs), 3, nan); 
            
            %keeping track of n_reps for current trial type
            n_reps_vec = [n_reps_vec; length(curr_trs)];
            
                            
%           olf2_od_n
%           plot(squeeze(PID_traces_matched(:, curr_trs, 1)))
            stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);  %computing stimulus on and off frame numbers for olf1 and olf2
            
            curr_mean_traces = squeeze(mean(dff_data_mat_f(:, union_sig_cells, curr_trs), 3, 'omitnan'));
            
            %computing response sizes and logging only for sig cells
            stim_frs_simp = stim_frs{2};
            stim_frs_simp(2) = stim_frs_simp(2) + round(2./frame_time);
            mean_resps = mean(mean(dff_data_mat(stim_frs_simp(1):stim_frs_simp(2), :, curr_trs), 3, 'omitnan'), 1);
            nonsig_cells = find(sig_cell_mat(:, stim_type_n) == 0);
            mean_resps(nonsig_cells) = 0;
            saved_resp_sizes(:, stim_type_n) = mean_resps';
            
            
            %Logging single trial responses to compute PC representations later on
            curr_resps_pca = squeeze(mean(dff_data_mat(stim_frs_simp(1):stim_frs_simp(2), :, curr_trs), 'omitnan'));
            curr_insig_cells = find(sig_cell_mat(:, stim_type_n) == 0);
            curr_resps_pca(curr_insig_cells, :) = 0;
            saved_resps_pca(:, :, stim_type_n) = curr_resps_pca;   
            
            %identifying mis-identified junk cells
            stim_on = stim_frs{1};
            stim_on = stim_on(1) - 1;
            bad_cells = mean(abs(curr_mean_traces(1:stim_on, :)), 1, 'omitnan');
            bad_cells = find(bad_cells > 0.50);
            bad_cells_all = [bad_cells_all, bad_cells];
            
            %logging current fly's response traces
            resp_trace_mat = pad_n_concatenate(resp_trace_mat, curr_mean_traces, 3, nan);            
        end
        
        %getting rid of misidentified junk cells
        resp_trace_mat(:, bad_cells_all, :) = [];
        bad_cells_all = [];
        
        %subtracting the appropriate single pulse responses from time-aligned transition response traces to get rid of pulse1
        %responses, if specified
        if subtract_pulse1_resps == 1
            if curr_stim_type_simp == 0         %case when simple stimuli were delivered on olf1
                
            elseif curr_stim_type_simp == 1     %case when simple stimuli were delivered on olf2
                %subtracting PA simple pulse resps
                PA_simp_row = find(sig_cell_mat_key(:, 3) == 1 & sig_cell_mat_key(:, 2) == 0);
                PA_pulse1_row = find(sig_cell_mat_key(:, 1) == 3 | sig_cell_mat_key(:, 1) == 4);      %PA odor number switched from 3 to 4 in 2021
                resp_trace_mat(:, :, PA_pulse1_row) = resp_trace_mat(:, :, PA_pulse1_row) - resp_trace_mat(:, :, PA_simp_row);
                
                %subtracting BA simple pulse resps
                BA_simp_row = find(sig_cell_mat_key(:, 3) == 3 & sig_cell_mat_key(:, 2) == 0);
                BA_pulse1_row = find(sig_cell_mat_key(:, 1) == 10 | sig_cell_mat_key(:, 1) == 9);      %BA odor number switched from 10 to 9 in 2021
                resp_trace_mat(:, :, BA_pulse1_row) = resp_trace_mat(:, :, BA_pulse1_row) - resp_trace_mat(:, :, BA_simp_row);
                
                
            else
            end
            
        else
        end
        
       
        
        %pooling sig cells across flies
        if pool_flies_traces == 1
            all_resp_traces = pad_n_concatenate(all_resp_traces, resp_trace_mat, 2, nan);
        elseif pool_flies_traces == 0 && dir_n == single_fly_n     %GKCs 5, A'\B' KCs 5, 
            all_resp_traces = pad_n_concatenate(all_resp_traces, resp_trace_mat, 2, nan);            
        else
        end
        
        saved_resp_sizes_all = [saved_resp_sizes_all; saved_resp_sizes];
        n_cells_all = [n_cells_all; n_cells];
        disp(curr_dir)
        
        %training logistic regression decoders
        %calling logistic regression main_function
        stim_frs_plt = stim_frs{1};
        try
            [decoder_resp_traces, pred_accuracy_mat] = do_log_regression(single_traces_mat, n_reps_vec, sig_cell_mat_key, stim_frs_plt);        
        catch
            keyboard
        end
        decoder_resp_traces = movmean(decoder_resp_traces, 5, 1);
        all_fly_accuracy_mat = [all_fly_accuracy_mat; mean(pred_accuracy_mat, 'omitnan')];
        %Plotting prob traces for single fly
        %Plotting mean decoder prediction probability traces as a function of time (all leave-one-out test predictions)
        color_vecs_decoders = [0, 0, 0; 0.65, 0.65, 0.65; [123,50,148]./256];
        for stim_type = 1:size(sig_cell_mat_key, 1)
            curr_resp_mat = decoder_resp_traces(:, :, stim_type);
            if stim_type <=3    %case for single pulse stimuli
                stim_type_n = stim_type;
                fig_n = 96 + stim_type_n;
                for trace_type = 1:3 
                    curr_od_traces = curr_resp_mat(:, ((trace_type - 1).*8 + 1):(trace_type.*8));
                    mean_prob_trace = mean(curr_od_traces, 2, 'omitnan');
                    se_prob_trace = std(curr_od_traces, [], 2, 'omitnan')./sqrt(size(curr_od_traces, 2));
                    figure(fig_n)                    
                    shadedErrorBar([], mean_prob_trace, se_prob_trace, {'Color', color_vecs_decoders(trace_type, :)});
                    hold on 
                end
                ylabel(['probability odor ' num2str(stim_type_n)]);
                set_xlabels_time(fig_n, frame_time, 10);
                fig_wrapup(fig_n, [], [25, 25], .6);
                add_stim_bar(fig_n, stim_frs_plt, color_vecs_decoders(stim_type_n, :));
                ax_vals = axis;
                plot([0, ax_vals(2)], [.5, .5], 'Color', [.8, .2, .2]);
                hold off
            elseif stim_type > 3    %case for transition stimuli
                stim_type_n = stim_type - 3;
                fig_n = 100 + stim_type_n;
                for trace_type = 1:2
                    trace_type_plt = [1, 2];
                    trace_type_plt(trace_type_plt == trace_type) = [];
                    curr_od_traces = curr_resp_mat(:, ((trace_type - 1).*8 + 1):(trace_type.*8));
                    mean_prob_trace = mean(curr_od_traces, 2, 'omitnan');
                    se_prob_trace = std(curr_od_traces, [], 2, 'omitnan')./sqrt(size(curr_od_traces, 2));
                    figure(fig_n)                    
                    shadedErrorBar([], mean_prob_trace, se_prob_trace, {'Color', color_vecs_decoders(trace_type_plt, :)});
                    hold on 
                end
                ylabel(['probability odor ' num2str(stim_type_n)]);
                set_xlabels_time(fig_n, frame_time, 10);
                fig_wrapup(fig_n, [], [25, 25], .6);
                stim_frs_plt = [stim_frs{1}; stim_frs{2}];
                if stim_type == 4
                    curr_color_vecs = [color_vecs_decoders(1, :); color_vecs_decoders(2, :)];
                elseif stim_type == 5
                    curr_color_vecs = [color_vecs_decoders(2, :); color_vecs_decoders(1, :)];
                else
                end
                add_stim_bar(fig_n, stim_frs_plt, curr_color_vecs);
                ax_vals = axis;
                plot([0, ax_vals(2)], [.5, .5], 'Color', [.8, .2, .2]);
                hold off                                                                                                                                                      
                
            else
            end
            
        end
               
        %computing and plotting PC representations for current fly
        ave_resp_mat = squeeze(mean(saved_resps_pca, 2, 'omitnan'));
        all_data_mat = reshape(saved_resps_pca, size(saved_resps_pca, 1), (size(saved_resps_pca, 2).*size(saved_resps_pca, 3)));
        PC_wts = pca(all_data_mat');
        
        PC1_resps = squeeze(sum(saved_resps_pca.*repmat(PC_wts(:, 1), 1, size(saved_resps_pca, 2), size(saved_resps_pca, 3)), 1, 'omitnan'));
        PC2_resps = squeeze(sum(saved_resps_pca.*repmat(PC_wts(:, 2), 1, size(saved_resps_pca, 2), size(saved_resps_pca, 3)), 1, 'omitnan'));
        
%         mean_PC1_resps = squeeze(sum(ave_resp_mat.*repmat(PC_wts(:, 1), 1, size(ave_resp_mat, 2)), 1, 'omitnan'));
%         mean_PC2_resps = squeeze(sum(ave_resp_mat.*repmat(PC_wts(:, 2), 1, size(ave_resp_mat, 2)), 1, 'omitnan'));
        
        color_vecs = [PA_color; BA_color; EL_color];
        %plotting in PC space
        for od_n = 1:3
            figure(14)
            plot(PC1_resps(:, od_n), PC2_resps(:, od_n), '.', 'Color', color_vecs(od_n, :), 'markerSize', 15.3);
            hold on
            %plot(mean_PC1_resps(od_n), PC2_resps(od_n), '.', 'Color', color_vecs(od_n, :).*0.75, 'markerSize', 24);
            
        end
        xlabel('PC1')
        ylabel('PC2')
        hold off
        fig_wrapup(14,[], [25, 25], 0.6)
        if pause_PCAs == 1 && dir_n == single_fly_n
            keyboard
        else
        end
        saved_resps_pca = [];
     
        %keyboard
    end
     
         
     od_name_vec = [];
     q_resp_mat = [];
     
     %making sure BA-PA trials are analyzed before PA-BA trials
     sig_cell_mat_key_orig = sig_cell_mat_key;
     sig_cell_mat_key(4, :) = sig_cell_mat_key_orig(5, :);
     sig_cell_mat_key(5, :) = sig_cell_mat_key_orig(4, :);
     
     for stim_type_n = 1:size(sig_cell_mat_key, 1)
        curr_stim_vec = sig_cell_mat_key(stim_type_n, :);

        olf1_od_n = curr_stim_vec(1, 1);
        olf1_dur = curr_stim_vec(1, 2);
        olf2_od_n = curr_stim_vec(1, 3);
        olf2_dur = curr_stim_vec(1, 4);

        if olf1_dur ~= 0 && olf2_dur ~= 0
            curr_stim_type = 2;     %handover stimuli
        elseif olf2_dur == 0
            curr_stim_type = 0;     %simple stimulus, delivered on olf1
        elseif olf1_dur == 0
            curr_stim_type = 1;     %simple stimulus, delivered on olf2                
        end
        
        
        %piggy-backing on last dataset to determine stim_frs, assuming they
        %are the same across all datasets for the current stimulus type
        curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == olf1_od_n & stim_mat_simple(:, dur_col_ns(1)) == olf1_dur &... 
                                stim_mat_simple(:, od_col_ns(2)) == olf2_od_n & stim_mat_simple(:, dur_col_ns(2)) == olf2_dur);
            
        stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);  %computing stimulus on and off frame numbers for olf1 and olf2
        
        %constructing odor stimulus name string
        olf1_od_ind = find(odor_list_olf1 == olf1_od_n);
        olf1_od_name = od_names{olf1_od_ind};
        olf2_od_ind = find(odor_list_olf2 == olf2_od_n);
        olf2_od_name = od_names{olf2_od_ind};

        if curr_stim_type == 0  %simple stimulus on olf1
            curr_stim_name = [olf1_od_name, ' - single pulse'];
            curr_save_name = [olf1_od_name, '_single_pulse_', KC_type];
            stim_frs = stim_frs{1};z                      
            curr_color = color_vec(olf1_od_ind, :);
        elseif curr_stim_type == 1 %simple stimulus on olf2
            curr_stim_name = [olf2_od_name, ' - single pulse'];
            curr_save_name = [olf2_od_name, '_single_pulse_', KC_type];
            stim_frs = stim_frs{2};
            curr_color = color_vec(olf2_od_ind, :);
        elseif curr_stim_type == 2 %odor transition stimulus
            curr_stim_name = [olf1_od_name, ' - ', olf2_od_name, ' transition'];
            curr_save_name = [olf1_od_name, olf2_od_name, '_transition_', KC_type];
            stim_frs = [stim_frs{1}; stim_frs{2}];
            curr_color = [color_vec(olf1_od_ind, :); color_vec(olf2_od_ind, :)];
        else
        end
        
        curr_mean_traces = squeeze(all_resp_traces(:, :, stim_type_n));
        
        %sorting cells by response size to first stim_type
        if stim_type_n == 1 || stim_type_n == 4
            wt_vec = mean(curr_mean_traces(stim_frs(1):stim_frs(2), :), 1, 'omitnan');     %response sizes on stim_type 1
        else
        end
        %sorting by wt_vec
        curr_mean_traces = [wt_vec; curr_mean_traces]';
        curr_mean_traces = sortrows(curr_mean_traces)';
        curr_mean_traces(1, :) = [];
        
        %computing mean resp over odor period (pulse2 for transition pulses)
        if curr_stim_type == 1
            q_resp_vec = max(curr_mean_traces(stim_frs(1, 1):(stim_frs(1, 2) + round(2./frame_time)), :));      %simple stimuli
            od_name_vec = [od_name_vec; {[olf2_od_name, ' simple']}];
        elseif curr_stim_type == 2
            q_resp_vec = max(curr_mean_traces(stim_frs(1, 2):(stim_frs(2, 2) + round(2./frame_time)), :));      %transition stimuli, averaging only pulse2 period
            od_name_vec = [od_name_vec; {[olf2_od_name, ' transition']}];
        else
        end
        
        q_resp_mat = [q_resp_mat; q_resp_vec];
       
        
        fig_h = figure('Name', curr_stim_name);
        
        if isempty(findstr(curr_dir_list_path, 'd5HT1b')) == 0
            high_val = 1.5;
        elseif isempty(findstr(curr_dir_list_path, 'c305a')) == 0
            high_val = 6;
        elseif isempty(findstr(curr_dir_list_path, 'c739')) == 0
            high_val = 3;
        else
        end
        imagesc(curr_mean_traces', [0, high_val]);
        colormap(greymap);
        set_xlabels_time(fig_h, frame_time, 15);
        ylabel('cell number');
        fig_wrapup_mod(fig_h, 'tall', []);
        fig_wrapup(fig_h, [], [25, 15], 0.6)
        add_stim_bar(fig_h, stim_frs, curr_color);
        
        write_data_cols = curr_mean_traces;
        if curr_stim_type < 2
            xls_path = [paper_save_dir,  curr_save_name, '_traces.xls'];
        elseif curr_stim_type == 2
            xls_path = [paper_save_dir_tr,  curr_save_name, '_traces.xls'];
        end
        [c] = write_xls_header([], write_data_cols, xls_path);
        write_data_cols = [];
       
    end
    
    %plotting scatter plots and computing corrcoefs
    [corr_mat, p_mat] = corrcoef(q_resp_mat');
    
    pair_list = [1, 2;...   %1. Simple PA - Simple BA
                 1, 3;...   %2. Simple PA - Simple EL
                 
                 2, 3;...   %3. Simple BA - Simple EL
                 5, 4;...   %4. Transition PA - Transition BA
                 1, 5;...   %4. Simple PA - Transition BA
                 2, 4;...   %4. Simple BA - Transition PA
                 1, 4;...   % Simple PA - Transition PA Note: This is only meaningful if pulse1 responses are subtracted from the transition trial responses
                 2, 5;...   % Simple BA - Transition BA Note: This is only meaningful if pulse1 responses are subtracted from the transition trial responses
                 ];
    
    for pair_n = 1:size(pair_list, 1)
        curr_pair = pair_list(pair_n, :);
        fig_h = figure('Name', [od_name_vec{curr_pair(1)}, ' - ', od_name_vec{curr_pair(2)}]);
        plot(q_resp_mat(curr_pair(1), :), q_resp_mat(curr_pair(2), :), 'O', 'MarkerSize', 4, 'MarkerEdgeColor', [0.65, 0.65, 0.65])
        xlabel([od_name_vec{curr_pair(1)}, ' responses']);
        ylabel([od_name_vec{curr_pair(2)}, ' responses']);
        ax_vals = axis;
        max_ax = max(ax_vals);
        axis([-0.2, max_ax, -0.2, max_ax]);
        ax_vals = axis;
        x_val = ax_vals(1) + ax_vals(2).*0.05;
        y_val = ax_vals(4).*0.95;
        r_val = num2str(round(corr_mat(curr_pair(1), curr_pair(2)), 3));
        p_val = num2str(round(p_mat(curr_pair(1), curr_pair(2)), 3));
        if p_val == '0'
            p_val = '< 0.001';
        else
        end
        text(x_val, y_val, ['r ', r_val, ', p ', p_val])
        fig_wrapup(fig_h, [], [25, 30], 0.6);
        
        %saving data underlying plots
        write_data_cols = [(q_resp_mat(curr_pair(1), :))', (q_resp_mat(curr_pair(2), :))'];
        header_row = [{od_name_vec{curr_pair(1)}}, {od_name_vec{curr_pair(2)}}];
        if pair_n == 1            
            xls_path = [paper_save_dir_sfig,  'scatter_cell_resps_similar', KC_type, '.xls'];
            
            %writing data behind plot to file
            [c] = write_xls_header(header_row, write_data_cols, xls_path);
            write_data_cols = [];
        elseif pair_n == 2
            xls_path = [paper_save_dir_sfig,  'scatter_cell_resps_dissimilar', KC_type, '.xls'];
            
            %writing data behind plot to file
            [c] = write_xls_header(header_row, write_data_cols, xls_path);
            write_data_cols = [];
            
        else
        end
        
        
    end
    
        
   
    %computing mean Euclidean distances between odors for each fly 
    dists_saved = zeros(size(n_cells_all, 1), 4) + nan;
    for fly_n = 1:size(n_cells_all, 1)
        if fly_n > 1
            curr_celln_st = sum(n_cells_all(1:(fly_n - 1))) + 1; 
        else
            curr_celln_st = 1;
        end
        curr_celln_end = sum(n_cells_all(1:fly_n));
        
        curr_resp_mat = saved_resp_sizes_all(curr_celln_st:curr_celln_end, :)';
        pad = zeros(1, size(curr_resp_mat, 2));
        curr_resp_mat_p = [pad; curr_resp_mat];   %adding a point for the origin        
        dist_mat = pdist(curr_resp_mat_p);
        dist_mat = squareform(dist_mat);
        %normalizing all distances to distance of EL from origin
        dist_mat = dist_mat./(dist_mat(1, 2));
        dists_saved(fly_n, :) = [dist_mat(2, 3), dist_mat(2, 4), dist_mat(3, 4), dist_mat(5, 6)];     %PA-BAsim, PA-ELsim, BA-ELsim, BAPA - PABA
        
    end
    
    
 
    
    %plotting Euclidean distances
    fig_h = figure();
    marker_colors = repmat([0.6, 0.6, 0.6], 4, 1);
    col_pairs = [];
    line_colors = marker_colors;
    mean_color = [0 ,0, 0];
    xlabels = [{'PA-BA', 'PA-EL', 'BA-EL', 'BAPA-PABA'}];
    fig_h = scattered_dot_plot_ttest(dists_saved(:, 1:3), fig_h, 2.5, 4, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean');
    ylabel('norm. Euclidean dist.')
    ax_vals = axis;
    ax_vals(4) = 4;
    axis(ax_vals);
    fig_wrapup(fig_h, [], [25, 25], 0.6);
    
    p_simdif1 = signrank(dists_saved(:, 1), dists_saved(:, 2))
    p_simdif2 = signrank(dists_saved(:, 1), dists_saved(:, 3))
    p_simdif1_2_corr = bonf_holm([p_simdif1, p_simdif2], 0.05)
    
    
    
    %plotting Euclidean distances with transition resp distances
    fig_h = figure();
    marker_colors = repmat([0.6, 0.6, 0.6], 4, 1);
    col_pairs = [];
    line_colors = marker_colors;
    mean_color = [0, 0, 0];
    xlabels = [{'PA-BA', 'BP-PB'}];
    fig_h = scattered_dot_plot_ttest(dists_saved(:, [1,4]), fig_h, 2.5, 4, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean');
    ylabel('norm. Euclidean dist.')
    ax_vals = axis;
    ax_vals(4) = 2;
    axis(ax_vals);
    fig_wrapup(fig_h, [], [25, 25], 0.6);
    
    p_simdiftr = signrank(dists_saved(:, 1), dists_saved(:, 4))
    
    
    %plotting logistic regression decoder accuracies for simple trials
    fig_h = figure()
    xlabels = [{'PA', 'BA', 'EL'}];
    fig_h = scattered_dot_plot_ttest(all_fly_accuracy_mat(:, 1:3).*100, fig_h, 2.5, 4, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean');
    ylabel('decoder accuracy')
    ax_vals = axis;
    ax_vals(4) = 100;
    axis(ax_vals);
    hold on
    plot([0, ax_vals(2)], [50, 50], 'Color', [1, 0, 0]);
    fig_wrapup(fig_h, [], [15, 30], 0.6);
    write_data_cols = all_fly_accuracy_mat(:, 1:3);
    xls_path = [paper_save_dir, KC_type, '_logistic_decoder_acc.xls'];
    col_heads = [{'PA'}, {'BA'}, {'EL'}];
    [c] = write_xls_header(col_heads, write_data_cols, xls_path);
    write_data_cols = [];
   
    p_decoder_s_PAEL = signrank(all_fly_accuracy_mat(:, 1), all_fly_accuracy_mat(:, 3))
    p_decoder_s_BAEL = signrank(all_fly_accuracy_mat(:, 2), all_fly_accuracy_mat(:, 3))
    p_decoder_s_corr = bonf_holm([p_decoder_s_PAEL, p_decoder_s_BAEL], 0.05)
    
    %plotting logistic regression decoder accuracies for simple and transition trials
    fig_h = figure()
    xlabels = [{'PA', 'BA-PA', 'BA', 'PA-BA'}];
    fig_h = scattered_dot_plot_ttest(all_fly_accuracy_mat(:, [1, 5, 2, 4]).*100, fig_h, 2.5, 4, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean');
    ylabel('decoder accuracy')
    ax_vals = axis;
    ax_vals(4) = 100;
    axis(ax_vals);
    hold on
    plot([0, ax_vals(2)], [50, 50], 'Color', [1, 0, 0]);
    fig_wrapup(fig_h, [], [25, 25], 0.6);
    
    p_decoder_trtr = signrank(all_fly_accuracy_mat(:, 4), all_fly_accuracy_mat(:, 5))
    p_decoder_simtrPA = signrank(all_fly_accuracy_mat(:, 1), all_fly_accuracy_mat(:, 5))
    p_decoder_simtrBA = signrank(all_fly_accuracy_mat(:, 2), all_fly_accuracy_mat(:, 4))
    p_decoder_tr_corr = bonf_holm([p_decoder_simtrPA, p_decoder_simtrBA], 0.05)
    
    
end



