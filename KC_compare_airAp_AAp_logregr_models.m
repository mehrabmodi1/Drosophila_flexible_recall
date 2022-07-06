clear all
close all

dataset_list_paths = [ ...
                         {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\d5HT1b_similar_od_handovers.xls'};...
                         {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\c305a_similar_od_handovers.xls'};...
                         {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\c739_similar_od_handovers.xls'};...
                      ];

[del, odor_names1] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
[del, odor_names2] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);
paper_save_dir = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\Fig5_KCs_transitions\';

odor_names2{3} = 'Butyl acetate';
od_names = [{'PA'}, {'BA'}, {'EL'}];

paired_color = [0,136,55]./256;
unpaired_color = [166,219,160]./256;
EL_color = [123,50,148]./256;
PA_color = [0,136,55]./256;
BA_color = [166,219,160]./256;
mean_color = [0.8, 0.4, 0.4];

global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;

re_extract_KC_data = 0;
plotscore_flies = 0;    %Switches between plotting and testing indiv models or model outputs averaged for each fly


%Fit related parameters
multi_fly_MBON = 1;     %This switches between fitting KC trials from a single fly to resampled MBON responses from multiple flies or the single trials from a single fly.
fit_simp_and_trans = 1; %This switches between (0) fitting log regressors to simple trials only (simulate biology) or (1) simple AND transition trials (most generous test of hypothesis).
KC_act_threshold = 0; %This floors all KC activity values below this threshold during fitting and decoding
training_odor_vec = [1, 2, 3, 4, 5];    %train on all odors (NOTE: Always tested on all odors)
%training_odor_vec = [1, 3, 4, 5];    %ignore CS- only    (NOTE: Always tested on all odors)
%training_odor_vec = [1, 2, 4, 5];    %ignore Dissimilar only    (NOTE: Always tested on all odors)
%training_odor_vec = [1, 2, 3];    %ignore transitions    (NOTE: Always tested on all odors)
%training_odor_vec = [1, 3];      %only train on paired and dissimilar odors.
%training_odor_vec = 1;           %only train on paired odor.   results in trivial training where all outputs are 0.
%training_odor_vec = [2, 5];


%initialize weights with an offset to affect ignored odor decoding
%wt_lookup_vec = [2, 1, 1.5];
wt_lookup_vec = [1, 1, 1];     %default
lambda = 0.0;         %regularization constant


%Use 2s response integration window
%KC_win = 2;     %deault - 7 use KC response data extracted over KC_win s after relevant odor pulse onset
KC_win = 7;


if KC_win == 7
    save_path_base = 'C:\Data\Data\Analysed_data\Analysis_results\KC_transition_logregr\';
elseif KC_win == 2
    save_path_base = 'C:\Data\Data\Analysed_data\Analysis_results\KC_transition_logregr_2s\';
end

fig_save_path = 'C:\Users\Mehrab\Dropbox (HHMI)\data_sharing\figure_for_Glenn_20220314\KC_analysis_figs\logistic_regr_models\';

%extracting, re-formatting and writing KC response data to disk
if re_extract_KC_data == 1
    for list_n = 1:size(dataset_list_paths, 1)
        curr_dir_list_path = dataset_list_paths{list_n, 1};
        curr_dir_list_path = update_list_path(curr_dir_list_path);
        [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
        n_dirs = size(dir_list, 1);
        dataset_list_name = findstr(curr_dir_list_path, 'list_');
        dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));

        %parsing path to identify current KC type
        slashi = findstr(curr_dir_list_path, '\');
        slashi = slashi(length(slashi));
        uscorei = findstr(curr_dir_list_path, '_');
        uscorei = uscorei - slashi;
        uscorei = uscorei(uscorei > 0);
        KC_type = curr_dir_list_path((slashi + 1):(slashi + uscorei - 1));


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
            [del, all_bad_trs] = cell_data_quality_control(dff_data_mat_f, stim_mat, stim_mat_simple, column_heads, sig_cell_mat, sig_cell_mat_key, 1, frame_time);

            union_sig_cell_vec = sum(sig_cell_mat, 2);
            union_sig_cell_vec = sign(union_sig_cell_vec);
            union_sig_cells = find(union_sig_cell_vec == 1);
            %sum(sig_cell_mat)./size(sig_cell_mat, 1);
            union_nonsig_cells = 1:size(sig_cell_mat, 1);
            union_nonsig_cells(union_sig_cells) = [];

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

                curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == olf1_od_n & stim_mat_simple(:, dur_col_ns(1)) == olf1_dur &... 
                                    stim_mat_simple(:, od_col_ns(2)) == olf2_od_n & stim_mat_simple(:, dur_col_ns(2)) == olf2_dur);

                %logging traces sorted by stim type to use for logistic regression decoding.
                single_traces_mat = pad_n_concatenate(single_traces_mat, dff_data_mat_f(:, union_sig_cells, curr_trs), 4, nan); 

                stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);  %computing stimulus on and off frame numbers for olf1 and olf2

                %computing response sizes and logging only for sig cells
                stim_frs_olf2 = stim_frs{2};
                stim_frs_olf2(2) = stim_frs_olf2(1) + round(KC_win./frame_time);
                
                %simple and transition stimuli, both delivered with olf2
                t_win_resps = squeeze(mean(dff_data_mat(stim_frs_olf2(1):stim_frs_olf2(2), union_sig_cells, curr_trs), 1, 'omitnan'));
                               
                saved_resp_sizes = pad_n_concatenate(saved_resp_sizes, t_win_resps, 3, nan);
                
            end

            fly_data.resp_traces = single_traces_mat;
            fly_data.resp_sizes = saved_resp_sizes;
            fly_data.stim_mat_key = sig_cell_mat_key;
            fly_data.stim_frs = stim_frs;
            fly_data.stim_mat = stim_mat_simple;
            fly_data.KC_type = KC_type;
            fly_data.sig_cell_mat = sig_cell_mat(union_sig_cells, :);
            
            save_dir = [save_path_base, KC_type, '\fly', num2str(dir_n)];
            mkdir(save_dir);
            save([save_dir, '\fly_data.mat'], 'fly_data');
            
        end

       



    end
else
end


%Reading in KC and MBON response data
%reading in MBON data
% MBON_path = [save_path_base, 'MBON\'];
% MBON_simp_data = load([MBON_path, 'simple_stim_data.mat']);
% MBON_simp_data = MBON_simp_data.all_MBON_data;
% MBON_simp_paired_ods = MBON_simp_data.paired_ods_olf1;
% MBON_simp_resps_all = MBON_simp_data.resp_mat;
% 
% MBON_tr_data = load([MBON_path, 'transition_stim_data.mat']);
% MBON_tr_data = MBON_tr_data.all_MBON_data;
% MBON_tr_paired_ods = MBON_tr_data.paired_ods_olf1;
% MBON_tr_resps_all = MBON_tr_data.resp_mat;

%reading in KC data, doing fit and logging results
KC_types = [{'d5HT1b'}, {'c305a'}, {'c739'}];

for KC_type_n = 1:3
    fig_n = 0;
    KC_type = KC_types{KC_type_n};
    n_flies = length(dir([save_path_base, KC_type])) - 2;
    
    
    %looking up initial weight matrix offset from manually specified vector of offsets for each KC subtype
    initial_wts_offset = wt_lookup_vec(KC_type_n);
        
    %result log variables
    test_acc_s_all = [];
    test_transition_acc_all = [];
    test_tr_vals_all = [];    
    wts_all = [];
    resp_vecs_all = [];
    median_wts = [];
    median_resps_all = [];
    test_tr_vals_fly = [];
    test_tr_vals_fly_ave = [];
    test_tr_vals_fly_pre = [];
    test_tr_vals_fly_ave_pre = [];
    for fly_n = 1:n_flies
        
        %skipping flies with few cells
        if KC_type_n <=2 && fly_n == 4  
            continue
        else
        end
        
        KC_data = load([save_path_base, KC_type, '\fly', num2str(fly_n), '\fly_data.mat']);
        KC_data = KC_data.fly_data;
        KC_resp_data = KC_data.resp_sizes;
        KC_resp_data(KC_resp_data < KC_act_threshold) = 0;
        
        fly_wts = [];
        fly_resps = [];
        %looping through PA and BA as the paired odors
        for od_ni = 1:2
            if od_ni == 1
                paired_od = 3;
            else
                paired_od = 1;
            end
            
            %picking up only simple or simple and transition KC resp data for training, as specified.
            if fit_simp_and_trans == 0
                KC_resp_data_use = KC_resp_data(:, :, 1:3);
            elseif fit_simp_and_trans == 1
                KC_resp_data_use = KC_resp_data(:, :, 1:5);
            else
            end
                
            %loop to train on all trials but 1 to use leave one out cross-validation
            for l_trial_n = 1:size(KC_resp_data, 2)
                curr_trs = 1:size(KC_resp_data, 2);
                curr_trs(l_trial_n) = [];
                inputs_train = KC_resp_data_use(:, curr_trs, :);
                inputs_test = KC_resp_data_use(:, l_trial_n, :);        %This test is run with simple pulse KC responses, to measure goodness of fit.

                outputs_train = zeros(length(curr_trs), 5);              %constructing perfect output vectors
                outputs_train(:, 3) = 1;                    %assigning EL responses as always high
                %outputs_train(:, 2) = 1;                    %temporarily assigning A' responses as always high
                outputs_test = zeros(1, 5);
                outputs_test(3) = 1;                    %setting EL as high

                %if fitting trans and simple responses, assigning each odor as high in turn
                %in Original KC data: 1 - Air-PA, 2 - Air-BA, 3 - Air-EL, 4 - PA-BA, 5 - BA-PA
                %Paired-od convention: 1 - Air-A, 2 - Air-A', 3 - Air-B, 4 - A'-A, 5 - A-A'. Hence have to re-arrange as below.
                inputs_train_orig = inputs_train;
                inputs_test_orig = inputs_test;
                if od_ni == 1   %PA is paired odor
                    inputs_train(:, :, 4) = inputs_train_orig(:, :, 5); %1 is Air-PA and 4 is now BA-PA 
                    inputs_train(:, :, 5) = inputs_train_orig(:, :, 4); %2 is Air-BA and 5 is now PA-BA
                    
                    inputs_test(:, :, 4) = inputs_test_orig(:, :, 5);
                    inputs_test(:, :, 5) = inputs_test_orig(:, :, 4);
                elseif od_ni == 2   %BA is paired odor
                    inputs_train(:, :, 1) = inputs_train_orig(:, :, 2);
                    inputs_train(:, :, 2) = inputs_train_orig(:, :, 1);
                    
                    inputs_test(:, :, 1) = inputs_test_orig(:, :, 2);
                    inputs_test(:, :, 2) = inputs_test_orig(:, :, 1);                    
                end
                
                outputs_train(:, 5) = 1;    %since transition odors were randomly re-arranged above, assigning a fixed position as high
                outputs_test(5) = 1;        

                %using only selected odors for training
                inputs_train  = inputs_train(:, :, training_odor_vec);
                outputs_train = outputs_train(:, training_odor_vec);
                
                
                %normalizing input vectors
                inputs_train(inputs_train > 5) = 5;
                max_x = max(max(inputs_train));
                inputs_train = inputs_train./max_x;

                inputs_test(inputs_test > 5) = 5;
                max_x = max(max(inputs_test));
                inputs_test = inputs_test./max_x;
                

                try
                    [wts, test_tr_vals, test_acc_simp, test_acc_transition, test_tr_vals_pre] = log_regr_MBON_model(inputs_train, outputs_train, inputs_test, outputs_test, initial_wts_offset, lambda);    %test here is run with left out simple pulse KC responses, to measure goodness of fit.     
                catch
                    wts = [];
                   
                end


                %logging results
                if isempty(wts) == 0
                    test_tr_vals_all = pad_n_concatenate(test_tr_vals_all, test_tr_vals, 2, nan);
                    test_acc_s_all = [test_acc_s_all; test_acc_simp];
                    test_transition_acc_all = [test_transition_acc_all; test_acc_transition];
                    fly_wts = pad_n_concatenate(fly_wts, wts, 2, nan);
                    median_response_vec = median(inputs_train, 2, 'omitnan');        %logging median across fitted KC response vectors 
                    fly_resps = pad_n_concatenate(resp_vecs_all, median_response_vec, 2, nan);
                    test_tr_vals_fly = pad_n_concatenate(test_tr_vals_fly, test_tr_vals, 2, nan);
                    test_tr_vals_fly_pre = pad_n_concatenate(test_tr_vals_fly_pre, test_tr_vals_pre, 2, nan);
                else
                    
                end
                
            end

        
        end
        if isempty(fly_wts) == 0
            wts_all = pad_n_concatenate(wts_all, fly_wts, 2, nan);
            median_wts = pad_n_concatenate(median_wts, median(fly_wts, 2), 2, nan);
            resp_vecs_all = pad_n_concatenate(resp_vecs_all, fly_resps, 2, nan);
            median_resps_all = pad_n_concatenate(median_resps_all, median(fly_resps, 2), 2, nan);
            
            test_tr_vals_fly_ave = pad_n_concatenate(test_tr_vals_fly_ave, mean(test_tr_vals_fly, 2, 'omitnan'), 2, nan);
            test_tr_vals_fly = [];
            test_tr_vals_fly_ave_pre = pad_n_concatenate(test_tr_vals_fly_ave_pre, mean(test_tr_vals_fly_pre, 2, 'omitnan'), 2, nan);
            test_tr_vals_fly_pre = [];
        else
        end
        
        
    end
    wts_all(1, :) = nan;
    median_wts(1, :) = nan;
        
    fig_h = figure('Name', 'Regressor output on test trials');
    xlabels = [{'A'}, {'A'''}, {'B'}, {'A''-A'}, {'A-A'''}];
    %color_vecs = repmat([0.7, 0.7, 0.7; 0, 0, 0], 5, 1);
    color_vecs = repmat([0, 0, 0], 5, 1);
    
    %interleaving all_wts_1 regressor outputs and fitted_wts regressor outputs for plot
    %test_tr_vals_fly_ave = test_tr_vals_fly_ave';
    test_tr_vals_fly_ave = test_tr_vals_fly_ave';
    %test_tr_vals_fly_ave_pre = test_tr_vals_fly_ave_pre';
    test_tr_vals_fly_ave_pre = [];
    test_tr_vals_all = test_tr_vals_all';
    
    plot_mat = [];
    for od_n = 1:size(test_tr_vals_fly_ave, 2)
        if plotscore_flies == 1        
            %plot_mat = [plot_mat, test_tr_vals_fly_ave_pre(:, od_n), test_tr_vals_fly_ave(:, od_n)];
            plot_mat = [plot_mat, test_tr_vals_fly_ave(:, od_n)];
        elseif plotscore_flies == 0
            %plot_mat = [plot_mat, test_tr_vals_fly_ave_pre(:, od_n), test_tr_vals_all(:, od_n)];
            plot_mat = [plot_mat, test_tr_vals_all(:, od_n)];
        else
        end
    end
    
    %plotting model performance
    fig_h = scattered_dot_plot_ttest(plot_mat, fig_h, .6, [3], 4, color_vecs, 0, [], [], xlabels, 2, [0.8, 0.5, 0.5], 2, 0.05, 0, 1, 'force_mean', [], 0, 4);
    hold on
    xlabel('decoded odor')
    ylabel('log. regr. output')
    colormap(greymap)
    fig_wrapup(fig_h, [], [50, 85], 0.6);
    ax_vals = axis;
    plot([ax_vals(1), ax_vals(2)], [0.5, 0.5], 'k--');
    
    write_data_cols = plot_mat;
    col_heads = xlabels;
    
    if length(training_odor_vec) == 3
        xls_path = [paper_save_dir,  'simple_trained_regressor_outputs_', KC_type, '.xls'];
    elseif length(training_odor_vec) == 5
        xls_path = [paper_save_dir,  'transition_trained_regressor_outputs_', KC_type, '.xls'];
    end
    [c] = write_xls_header(col_heads, write_data_cols, xls_path);
    write_data_cols = [];
    
    
%     fig_n = fig_n + 1;
%     fig_name = fig_h.Name;
%     savefig([fig_save_path, KC_type, ' ', num2str(fig_n), ' ', fig_name, '.fig']);
    
    
    %Statistical testing
    [p_simp, h] = signrank(plot_mat(:, 1), plot_mat(:, 2));
    [p_trans, h] = signrank(plot_mat(:, 4), plot_mat(:, 5));
    
    pvals = bonf_holm([p_simp, p_trans], 0.01)
    %pval = p_trans
    
    
    fig_h = figure('Name', 'weight vectors');
    imagesc(wts_all)
    colormap(greymap)
    ylabel('cell number')
    xlabel('test trial number (across flies)');
    title(KC_type)
    fig_n = fig_n + 1;
    fig_name = fig_h.Name;
    savefig([fig_save_path, KC_type, ' ', num2str(fig_n), ' ', fig_name, '.fig']);
    
       
    corr_vals_all = [];
    od_names = [{'A'}, {'A'''}, {'B'}, {'A''-A'}, {'A-A'''}];
    for train_od_ni = 1:length(training_odor_vec)
        train_od_n = training_odor_vec(train_od_ni);
        od_name = od_names{train_od_n};
        %relating weights to activity vector for each training odor
        fig_h = figure('Name', ['weights vs mean ', od_name, ' resps']);
        weights = reshape(median_wts, 1, []);
        A_activ = reshape(squeeze(median_resps_all(:, :, train_od_ni)), 1, []);
        A_activ = A_activ./max(A_activ);
        %plot(A_activ, ones(length(weights), 1).*initial_wts_offset, '.', 'Color', [0.7, 0.7, 0.7], 'MarkerSize', 6);
        %hold on
        plot(A_activ, weights, '.', 'Color', [0, 0, 0], 'MarkerSize', 6);
                
        [corr_vals, corr_p] = corrcoef(weights, A_activ, 'Rows', 'complete');  %computing correlation coeff
        corr_vals = corr_vals(1, 2);
        corr_p = corr_p(1, 2);
        ax_vals = axis();
        if corr_p > 0.001
            corr_label = ['R ', num2str(corr_vals), ' p ', num2str(round(corr_p, 3))];
        else
            corr_label = ['R ', num2str(corr_vals), ' p < 0.001'];
        end
        text((ax_vals(2).*0.2), (ax_vals(4).*0.8), corr_label);

        xlabel(['norm. resp. to ', od_name])
        ylabel('KC weight')
        fig_wrapup(fig_h, [], [5, 7], 0.6);
        corr_vals_all = [corr_vals_all; [corr_vals, corr_p]];
        fig_n = fig_n + 1;
        fig_name = fig_h.Name;
        savefig([fig_save_path, KC_type, ' ', num2str(fig_n), ' ', fig_name, '.fig']);
        
        %saving data underlying paper figures
        if train_od_n == 1      %only used weights v/s A in paper
            write_data_cols = [weights', A_activ'];
            
            %getting rid of nans
            del = isnan(write_data_cols(:, 1));
            write_data_cols(del, :) = [];
            del = isnan(write_data_cols(:, 2));
            write_data_cols(del, :) = [];
            
            col_heads = [{'weights'}, {'A_resps'}];

            if length(training_odor_vec) == 3
                xls_path = [paper_save_dir,  'simple_trained_A_resps_vs_wts_', KC_type, '.xls'];
            elseif length(training_odor_vec) == 5
                xls_path = [paper_save_dir,  'transition_trained_A_resps_vs_wts_', KC_type, '.xls'];
            end
            [c] = write_xls_header(col_heads, write_data_cols, xls_path);
            write_data_cols = [];
        else
        end
        
   end
   corr_vals_all
   
        
   keyboard
   close all
end




