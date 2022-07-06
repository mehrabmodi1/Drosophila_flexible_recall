clear all
close all

dataset_list_paths = [...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_simp_pairing_Berry.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_handover_pairing_Berry.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\MBONG2_PaBaEl_handover_pairing_Berry_noLEDctrl.xls'};...
                     ];
            
suppress_plots = 1;
plotting_quant_no_filt = 0;     %1 - only unfiltered traces used for all analysis and plotting - traces included. 0 - filtered traces used for everything.
cell_n = 1;

plot_diff_traces = 0;   %doesn't plot dF/F response traces, but post - pre response traces

global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);

script_name = mfilename;

y_ax_traces = 0.8;
y_ax_fit_traces = 0.6;

an_save_path = 'C:\Data\Data\Analysed_data\Analysis_results\PaBaEl_Gamma2\';
logregr_sv_path = 'C:\Data\Data\Analysed_data\Analysis_results\KC_transition_logregr\MBON\';
paper_save_dir = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\Fig4_MBON_transitions_fig\';
paper_save_dir_simple = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\Fig3_MBONG2Ap1_simple_fig\';
paper_save_dir_sfig = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\SFig3_4_MBON_transitions_fig\';
%to be used for no LED control data 
%paper_save_dir = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\SFig3_4_MBON_transitions_fig\noLED_control\';

force_resave = 1;

n_sec = 2;      %width of moving integration window in s

integ_win_s = 5;        %width of pulse integration window for response quantification
fly_n = 0;
for list_n = 1:size(dataset_list_paths, 1)
    paired_od_vec = [];
    
    curr_dir_list_path = dataset_list_paths{list_n, 1};
    curr_dir_list_path = update_list_path(curr_dir_list_path);
    [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
    n_dirs = size(dir_list, 1);
    dataset_list_name = findstr(curr_dir_list_path, 'list_');
    dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));
    
    if isempty(findstr(curr_dir_list_path, 'handover')) == 0
        if isempty(findstr(curr_dir_list_path, 'ELfirst')) == 0
            set_type_name = 'transitionsELfirst_MBONG2Ap1';
            set_list_type = 2;  %handover stimulus with EL on pulse1 case
        elseif isempty(findstr(curr_dir_list_path, 'ELsecond')) == 0
            set_type_name = 'transitionsELsecond_MBONG2Ap1';
            set_list_type = 3;  %handover stimulus with EL on pulse2 case
        elseif isempty(findstr(curr_dir_list_path, 'ELfirst')) == 1 && isempty(findstr(curr_dir_list_path, 'ELsecond')) == 1 && isempty(findstr(curr_dir_list_path, 'generaliz')) == 1
            set_type_name = 'transitions_MBONG2Ap1';
            set_list_type = 1;  %handover stimulus with no EL pulses
        elseif isempty(findstr(curr_dir_list_path, 'generaliz')) == 0
            set_list_type = 4;  %handover, generalization experiment stimuli (CS- with EL first or second)
            set_type_name = 'transitionsgeneraliz_MBONG2Ap1';
        else
        end
    else 
        set_list_type = 0;  %simple stimulus case
        set_type_name = 'simple_pulses_MBONG2Ap1';
    end
    
        
    saved_traces_all = [];
    saved_PID_traces_all = [];
    ctrast_resps_all = [];
    ctrst_t_win = 1;
    
    %loop to go through all experiment datasets listed in list file
    for dir_n = 1:n_dirs
        curr_dir = [dir_list{dir_n, 1}, '\'];
        curr_dir = manage_base_paths(curr_dir, 2)
       
        tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
        tif_times = tif_times.time_stamps;
        
        
        cd(curr_dir);
        tif_name = dir('*.tif');
        if isempty(tif_name) == 1
            continue
        else
        end
        fly_n = fly_n + 1;
        stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
        [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);
        [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig, PID_traces] = load_params_trains_modular(curr_dir, tif_times, frame_time);    %reading in trial stimulus parameters after matching time stamps to F traces
        
        odor_names1 = stim_mat.odourNames;
        odor_names2 = stim_mat.odourNames_olf2;
        PA_odn_olf1 = od_name_lookup(odor_names1, 'Pentyl acetate');
        PA_odn_olf2 = od_name_lookup(odor_names2, 'Pentyl acetate');
        BA_odn_olf1 = od_name_lookup(odor_names1, 'Butyl acetate');
        BA_odn_olf2 = od_name_lookup(odor_names2, 'Butyl acetate');
        EL_odn_olf1 = od_name_lookup(odor_names1, 'Ethyl lactate');
        EL_odn_olf2 = od_name_lookup(odor_names2, 'Ethyl lactate');
        
        if dir_n == 1
            ctrst_t_win = round(ctrst_t_win./frame_time);
        else
        end
        
        paired_color = [0,136,55]./256;
        unpaired_color = [166,219,160]./256;
        EL_color = [123,50,148]./256;
        mean_color = [0, 0, 0];
        %mean_color = [0.8, 0.4, 0.4];
        
        stim_mat_simple_nonans = stim_mat_simple;
        stim_mat_simple_nonans(isnan(stim_mat_simple)) = 0;
                
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
                
        %Reading in experimental parameters
        odor_list_olf1 = unique(stim_mat_simple(:, od_col_ns(1) ) );
        n_odors_olf1 = length(odor_list_olf1);
        odor_dur_list_olf1 = unique(stim_mat_simple(:, dur_col_ns(1) ) );
        n_od_durs_olf1 = length(odor_dur_list_olf1);
                
        odor_list_olf2 = unique(stim_mat_simple(:, od_col_ns(2) ) );
        n_odors_olf2 = length(odor_list_olf2);
        odor_dur_list_olf2 = unique(stim_mat_simple(:, dur_col_ns(2) ) );
        n_od_durs_olf2 = length(odor_dur_list_olf2);
              
        
        integ_win = round(integ_win_s./frame_time); %integration window for response quantification in frames
       
        %loading extracted raw fluorescence data matrices written by raw_dff_extractor
        raw_data_mat = load([curr_dir 'extracted_raw_data_mat.mat']);
        raw_data_mat = raw_data_mat.raw_data_mat;           %raw F traces extracted from ROIs
        raw_data_mat_orig = raw_data_mat;
        tif_n_col_n = find_stim_mat_simple_col('matched_tif_n', column_heads);
        n_cells = size(raw_data_mat, 2);
        if cell_n > n_cells
            cell_n = 1;
            disp('WARNING! specified cell_n not valid, reassigning cell_n = 1.');
        else
        end
        
        %inserting dummy trials (all nans) into raw_data_mat for pairing trials for
        %which no corress .tifs were acquired
        [raw_data_mat, good_tr_list] = match_up_rawmat_matchedtrs(raw_data_mat, stim_mat_simple, tif_n_col_n, good_tr_list);
                
        n_cells = size(raw_data_mat, 2);
        
        %calculating dF/F traces from raw data
        filt_time = 0.2;            %in s, the time window for boxcar filter for generating filtered traces
        [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir, good_tr_list);
        if plotting_quant_no_filt == 1
            dff_data_mat_f = dff_data_mat;
        else
        end
          

        dff_data_mat = dff_data_mat(:, cell_n, :);
        dff_data_mat_f  = dff_data_mat_f(:, cell_n, :);
        del = find(dff_data_mat_f < -1);
        dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones
                
       
        %identifying relevant odor numbers for each olfactometer
        pairing_tr_n = find(stim_mat_simple(:, led_on_col_n) == 1);
        
        if isempty(pairing_tr_n) == 1
            ctrl_set = 1;       %keeping track if this is a no LED ctrl dataset
            r_num = rand(1, 1);
            if r_num > 0.5
                %case for no LED ctrl datasets; paired odor assigned at random
                paired_od_n_olf2 = PA_odn_olf2;
                paired_od_n_olf1 = PA_odn_olf1;
            elseif r_num <= 0.5
                paired_od_n_olf2 = BA_odn_olf2;
                paired_od_n_olf1 = BA_odn_olf1;
            else
            end
        else
            ctrl_set = 0;
            
            prd_olf1_dur = stim_mat_simple(pairing_tr_n(1), dur_col_ns(1));
            if set_list_type == 0 || prd_olf1_dur > 1
                paired_od_n_olf1 = unique(stim_mat_simple(pairing_tr_n, od_olf1_col_n));
                if paired_od_n_olf1 == PA_odn_olf1
                    paired_od_n_olf2 = PA_odn_olf2;
                elseif paired_od_n_olf1 == BA_odn_olf1
                    paired_od_n_olf2 = BA_odn_olf2;
                else
                end
            elseif set_list_type > 0 && set_list_type < 4
                paired_od_n_olf2 = unique(stim_mat_simple(pairing_tr_n, od_olf2_col_n));
                
                if paired_od_n_olf2 == PA_odn_olf2
                    paired_od_n_olf1 = PA_odn_olf1;
                elseif paired_od_n_olf2 == BA_odn_olf2
                    paired_od_n_olf1 = BA_odn_olf1;
                else
                end
                
            elseif set_list_type == 4 %generalization dataset
                paired_od_n_olf2_orig = unique(stim_mat_simple(pairing_tr_n, od_olf2_col_n));
                if paired_od_n_olf2_orig == PA_odn_olf2
                    paired_od_n_olf2 = BA_odn_olf2;     %this is actually the CS- odor, assigned here to make subsequent parts of the script work
                    paired_od_n_olf1 = BA_odn_olf1;
                elseif paired_od_n_olf2_orig == BA_odn_olf2
                    paired_od_n_olf2 = PA_odn_olf2;     %this is actually the CS- odor, assigned here to make subsequent parts of the script work
                    paired_od_n_olf1 = PA_odn_olf1;     %this is actually the CS- odor, assigned here to make subsequent parts of the script work
                else
                end
            else
            end
        end
        
        paired_od_name = odor_names1{paired_od_n_olf1};
      
        if set_list_type ~= 4 
            if paired_od_n_olf1 == PA_odn_olf1
                unpaired_od_n_olf1 = BA_odn_olf1;
                unpaired_od_n_olf2 = BA_odn_olf2;
            elseif paired_od_n_olf1 == BA_odn_olf1
                unpaired_od_n_olf1 = PA_odn_olf1;
                unpaired_od_n_olf2 = PA_odn_olf2;
            else
            end
        elseif set_list_type == 4
           unpaired_od_n_olf1 = EL_odn_olf1;
           unpaired_od_n_olf2 = EL_odn_olf2;
        else
        end
       
        y_ax_lim = [];
        plot_means = 1;
        
        %manually specifying sets of pre and post pairing trial numbers
        if set_list_type <= 1
            tr_lists = [1:3; 7:9; 13:15];   %cases where EL alone trials exist
        elseif set_list_type > 1 && set_list_type ~= 4
            if stim_mat(9).trigger_scan == 0
                tr_lists = [1:2; 6:7; 11:12];   %cases where EL alone trials don't exist
            elseif stim_mat(9).trigger_scan == 1                   
                tr_lists = [1:2; 5:6; 10:11];   %cases where EL alone trials exist
            else
            end
        elseif set_list_type == 4   %generalization dataset
            tr_lists = [1:3; 7:9; 13:15];   %cases where EL alone trials exist
        else
        end
        
        if ctrl_set == 1
            tr_lists = [1:3; 4:6; 7:9];        %no LED stim, so no dummy, pairing trials inserted into stim_mat_simple
        else
        end
       
        
        %looping through pre, post1 and post2 trial types
        resp_mat = zeros(3, 3);
        
        saved_traces_curr = [];
        for tr_type_n = 1:2
            tr_list = tr_lists(tr_type_n, :);
            col_multiplier = 0.6.^(tr_type_n - 1);
            
            %paired odor response
            if set_list_type == 0
                curr_tr = find(stim_mat_simple(tr_list, od_olf1_col_n) == paired_od_n_olf1 & stim_mat_simple(tr_list, led_on_col_n) == 0);
                stim_frs_paired = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
                stim_frs = stim_frs_paired{1};      %integrating over entire pulse1 odor period for simple stim protocols
                stim_frs_pulse1 = stim_frs;
                stim_frs_EL = stim_frs_paired{1};
            elseif set_list_type == 1 || set_list_type == 2 || set_list_type == 4    %for handover without EL or handover EL first trials (also generalization trials with EL first)
                curr_tr = find(stim_mat_simple(tr_list, od_olf2_col_n) == paired_od_n_olf2 & stim_mat_simple(tr_list, led_on_col_n) == 0);
                stim_frs_paired = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
                stim_frs_pulse1 = stim_frs_paired{1};
                stim_frs = stim_frs_paired{2};      %integrating over entire pulse2 odor period for simple stim protocols
                stim_frs_EL = stim_frs_paired{1};
            elseif set_list_type == 3       %for handover, EL second trials
                curr_tr = find(stim_mat_simple(tr_list, od_olf1_col_n) == paired_od_n_olf1 & stim_mat_simple(tr_list, led_on_col_n) == 0);
                stim_frs_paired = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
                stim_frs_pulse1 = stim_frs_paired{1};
                stim_frs = stim_frs_paired{2};      %integrating over entire pulse2 odor period for simple stim protocols
                stim_frs_EL = stim_frs_paired{1};
            else
            end
            
            stim_frs(2) = stim_frs(2);  %adding on 2s after odor off to extend integration window
            stim_frs_EL(2) = stim_frs_EL(2); %adding on 2s after odor off to extend integration window
            curr_tr = tr_list(curr_tr);
            
            if length(curr_tr) ~= 1
                disp('more than one trial identified - please disambiguate (paired)')
                disp(['curr_tr = ' num2str(curr_tr)])
                figure(50)
                imagesc(stim_mat_simple, [0, 12]);
                curr_tr = input('Enter correct curr_tr: ');
            else
            end
            size(stim_mat_simple)
        
            curr_trace = squeeze(dff_data_mat_f(:, :, curr_tr));
            resp_vec(1, 1) = mean(mean(curr_trace(stim_frs(1):((stim_frs(1) + integ_win) + round(3./frame_time)) )));
            resp_vec_pulse1(1, 1) = mean(min(curr_trace((stim_frs_pulse1(2) - ctrst_t_win):((stim_frs_pulse1(2))) )));
            ctrst_vec_pulse2(1, 1) = mean(max(curr_trace(stim_frs(1):((stim_frs(1) + ctrst_t_win)) )));
            
            try
                saved_traces_curr(:, 1, tr_type_n) = curr_trace;
            catch
                keyboard
            end
                
            try
                saved_PID_traces_curr(:, 1, tr_type_n) = PID_traces(:, curr_tr);
            catch
                keyboard
            end
                
            fig_h = figure(1);
            set(fig_h, 'Position', [100, 100, 1200, 350]);
            subplot(1, 3, 1);
            plot(curr_trace, 'lineWidth', 1.5, 'Color', paired_color.*col_multiplier);
            hold on

            %un-paired odor response
            if set_list_type == 0
                curr_tr = find(stim_mat_simple(tr_list, od_olf1_col_n) == unpaired_od_n_olf1);
            elseif set_list_type == 1 || set_list_type == 2     %for handover without EL or handover EL first trials
                curr_tr = find(stim_mat_simple(tr_list, od_olf2_col_n) == unpaired_od_n_olf2);
            elseif set_list_type == 3       %for handover, EL second trials
                curr_tr = find(stim_mat_simple(tr_list, od_olf1_col_n) == unpaired_od_n_olf1);
            elseif set_list_type == 4       %for generalization, EL second trials
                curr_tr = find(stim_mat_simple(tr_list, od_olf2_col_n) == EL_odn_olf2);
            end
            
            curr_tr = tr_list(curr_tr);
            if length(curr_tr) ~= 1
                disp('more than one trial identified - please disambiguate (unpaired)')
                keyboard
                disp(['curr_tr = ' num2str(curr_tr)])
                figure(50)
                imagesc(stim_mat_simple, [0, 12]);
                curr_tr = input('Enter correct curr_tr: ');
            else
            end
            curr_trace = squeeze(dff_data_mat_f(:, :, curr_tr));
            resp_vec(1, 2) = mean(mean(curr_trace(stim_frs(1):((stim_frs(1) + integ_win) + round(3./frame_time)) )));
            resp_vec_pulse1(1, 2) = mean(min(curr_trace((stim_frs_pulse1(2) - ctrst_t_win):((stim_frs_pulse1(2))) )));
            ctrst_vec_pulse2(1, 2) = mean(max(curr_trace(stim_frs(1):((stim_frs(1) + ctrst_t_win)) )));
            try
                saved_traces_curr(:, 2, tr_type_n) = curr_trace;
            catch
                keyboard
            end
            
            saved_PID_traces_curr(:, 2, tr_type_n) = PID_traces(:, curr_tr);
            subplot(1, 3, 2)
            plot(curr_trace, 'lineWidth', 1.5, 'Color', unpaired_color.*col_multiplier);
            hold on

            %EL odor response or in the case of generalization expt, paired odor
            if set_list_type <= 1
                curr_tr = find(stim_mat_simple(tr_list, od_olf1_col_n) == 11);
                if isempty(curr_tr) == 1
                    curr_tr = find(stim_mat_simple(tr_list, od_olf2_col_n) == 4);
                else
                end
                curr_tr = tr_list(curr_tr);
                
                if length(curr_tr) ~= 1
                    disp('more than one trial identified - please disambiguate (EL)')
                    disp(['curr_tr = ' num2str(curr_tr)])
                    figure(50)
                    imagesc(stim_mat_simple, [0, 12]);
                    curr_tr = input('Enter correct curr_tr: ');
                else
                end
                
                curr_trace = squeeze(dff_data_mat_f(:, :, curr_tr));
                resp_vec(1, 3) = mean(mean(curr_trace(stim_frs_EL(1):((stim_frs_EL(1) + integ_win) + round(3./frame_time)) )));
                
            elseif set_list_type == 4   %generalization experiment trials, the single pulse trial here is a paired-odor trial
                curr_tr = find(stim_mat_simple(tr_list, od_olf2_col_n) == paired_od_n_olf2_orig);
                curr_tr = tr_list(curr_tr);
                curr_trace = squeeze(dff_data_mat_f(:, :, curr_tr));
                resp_vec(1, 3) = mean(mean(curr_trace(stim_frs_EL(1):((stim_frs_EL(1) + integ_win) + round(3./frame_time)) )));                
            else
                curr_trace = zeros(size(dff_data_mat_f, 1), 1) + nan;    %since EL alone trials don't exist, padding
                resp_vec(1, 3) = nan;
                
            end
            
            
            if set_list_type > 0
                stim_frs_bar_EL = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
                stim_frs_bar_EL = stim_frs_bar_EL{2};
            else
                stim_frs_bar_EL = compute_stim_frs_modular(stim_mat, curr_tr(1), frame_time);
                stim_frs_bar_EL = stim_frs_bar_EL{1};       %for simple stimulus dataset
            end
              
            saved_traces_curr(:, 3, tr_type_n) = curr_trace;
            saved_PID_traces_curr(:, 3, tr_type_n) = PID_traces(:, curr_tr);
            subplot(1, 3, 3)
            plot(curr_trace, 'lineWidth', 1.5, 'Color', EL_color.*col_multiplier);
            hold on
            resp_mat(tr_type_n, :) = resp_vec;
            resp_mat_pulse1(tr_type_n, :) = resp_vec_pulse1;
            ctrst_mat_pulse2(tr_type_n, :) = ctrst_vec_pulse2;
        end
        
        saved_traces_all = pad_n_concatenate(saved_traces_all, saved_traces_curr, 4, nan);
        saved_PID_traces_all = pad_n_concatenate(saved_PID_traces_all, saved_PID_traces_curr, 4, nan);
        clear saved_traces_curr
        clear saved_PID_traces_curr
        resp_mat_all(:, :, dir_n) = resp_mat;
        resp_mat_all_pulse1(:, :, dir_n) = resp_mat_pulse1;
        ctrst_mat_all_pulse2(:, :, dir_n) = ctrst_mat_pulse2;
        
        if suppress_plots == 0
            keyboard
        else
        end
        
        close figure 1
        
        paired_od_vec = [paired_od_vec; paired_od_n_olf2];
    end
    
    del = find(resp_mat_all > 10);
    resp_mat_all(del) = nan;
    del = find(saved_traces_all > 10);
    saved_traces_all(del) = nan;
    
    
    %saving MBON response data to file for separate analyses
    all_MBON_data.resp_mat = resp_mat_all(1:2, :, :);
    all_MBON_data.trace_mat = saved_traces_all;
    all_MBON_data.paired_ods_olf1 = paired_od_vec;
    
    if set_list_type == 0
        save_path = [logregr_sv_path, '\simple_stim_data.mat'];
    elseif set_list_type == 1
        save_path = [logregr_sv_path, '\transition_stim_data.mat'];
    else
    end
    save(save_path, 'all_MBON_data');
    
    
    
    %plotting averaged response traces
    %paired odor
    color_vecs = [0.65, 0.65, 0.65; 0, 0, 0; 0.65, 0.65, 0.65];
    figure(1)
    
    if set_list_type > 0
        stim_frs_bar = [stim_frs_paired{1}; stim_frs_paired{2}];
        %stim_frs_bar_EL = stim_frs_paired{2};
    elseif set_list_type == 0
        stim_frs_bar = stim_frs_paired{1};
        %stim_frs_bar_EL = stim_frs_paired{2};
    end
    
    table_labels = [{'pre'}; {'post'}];
    write_data_cols = [];
    for tr_type = 1:2
        if plot_diff_traces == 1
            diff_traces = saved_traces_all(1:(end - 5), 1, 2, :) - saved_traces_all(1:(end - 5), 1, 1, :);      %post traces - pre traces
            mean_trace = squeeze(mean(diff_traces, 4, 'omitnan'));
            se_trace = squeeze(std(diff_traces, [], 4, 'omitnan')./sqrt(size(saved_traces_all, 3)));
        else
            mean_trace = squeeze(mean(saved_traces_all(1:(end - 5), 1, tr_type, :), 4, 'omitnan')); 
            se_trace = squeeze(std(saved_traces_all(1:(end - 5), 1, tr_type, :), [], 4, 'omitnan')./sqrt(size(saved_traces_all, 3)));
        end
        plt_h = shadedErrorBar([], mean_trace, se_trace, {'Color', color_vecs(tr_type, :)}, 1);
        set(plt_h.edge(:), 'Color', 'none')
        if plot_diff_traces == 1
            break
        else
            hold on
        end
        
        write_data_cols = [write_data_cols, mean_trace, se_trace];  %logging data to write to file
                
    end
    hold off
    if plot_diff_traces == 1
        ylabel('paired odor responses (post - pre)')
    else
        ylabel('paired odor responses (dF/F)')
    end
    set_xlabels_time(1, frame_time, 10)
    ax_vals = axis;
    ax_vals(2) = 290;
    ax_vals(4) = 3;
    ax_vals(3) = -0.5;
    axis(ax_vals);
    fig_wrapup(1, [], [75, 90], 0.6)
    if set_list_type == 0
        add_stim_bar(1, stim_frs_bar, paired_color);
    elseif set_list_type == 1
        add_stim_bar(1, stim_frs_bar, [unpaired_color; paired_color]);
    elseif set_list_type == 2
        add_stim_bar(1, stim_frs_bar, [EL_color; paired_color]);
    elseif set_list_type == 3
        add_stim_bar(1, stim_frs_bar, [paired_color; EL_color]);
    elseif set_list_type == 4
        add_stim_bar(1, stim_frs_bar, [EL_color; unpaired_color]);      %in gen. datasets, the unpaired odor is coded as the paired odor in this analysis script
    else
    end
    
    if set_list_type == 0
        header_row = [{'A_pre_mean'}, {'A_pre_se'}, {'A_post_mean'}, {'A_post_se'}];
        xls_path = [paper_save_dir_simple,  'A_traces_', set_type_name, '.xls'];
    elseif set_list_type == 1
        header_row = [{'ApA_pre_mean'}, {'ApA_pre_se'}, {'ApA_post_mean'}, {'ApA_post_se'}];
        xls_path = [paper_save_dir,  'ApA_traces_', set_type_name, '.xls'];
    elseif set_list_type == 2
        header_row = [{'BA_pre_mean'}, {'BA_pre_se'}, {'BA_post_mean'}, {'BA_post_se'}];
        xls_path = [paper_save_dir,  'BA_traces_', set_type_name, '.xls'];
    elseif set_list_type == 3
        header_row = [{'AB_pre_mean'}, {'AB_pre_se'}, {'AB_post_mean'}, {'AB_post_se'}];
        xls_path = [paper_save_dir,  'AB_traces_', set_type_name, '.xls'];
    else
    end
    
    %writing data behind plot to file
    [c] = write_xls_header(header_row, write_data_cols, xls_path);
    write_data_cols = [];
    
    
    
    %fitting an exponential decay to mean trace
    if isempty(findstr(curr_dir_list_path, 'longpulse2')) ~= 1
        [del pki] = max(mean_trace(stim_frs(1):(stim_frs(2) + round(5./frame_time)) ));
        pki = pki + stim_frs(1);
        end_pt = min([(stim_frs(2) + round(20./frame_time)), size(mean_trace, 1)]);
        trace_sample = mean_trace(pki:end_pt);
        [x_vec, fit_params] = fit_exp_simple(trace_sample, frame_time);

        modelled_trace = exp(fit_params(2)).*exp(fit_params(1).*x_vec);
        disp(['tau = ',  num2str(-1.*(1./fit_params(1))), ' s']);
        
        figure(16)
        plt_h = shadedErrorBar([], mean_trace, se_trace, {'Color', color_vecs(tr_type, :)}, 1);
        set(plt_h.edge(:), 'Color', 'none')
        hold on
        x_vec = pki:1:(pki + size(modelled_trace, 1) - 1)';
        plot(x_vec, modelled_trace, ':r', 'lineWidth', 2);
        ax_vals = [0, 290, -1, 2];
        axis(ax_vals);
        ylabel('paired odor responses (dF/F)')
        set_xlabels_time(15, frame_time, 10)
        fig_wrapup(16, [], [75, 90], 0.6)
        add_stim_bar(16, stim_frs_bar, [unpaired_color; paired_color]);
    else
    end
    
    
    %unpaired odor
    figure(2)
    for tr_type = 1:2
        if plot_diff_traces == 1
            diff_traces = saved_traces_all(1:(end - 5), 2, 2, :) - saved_traces_all(1:(end - 5), 2, 1, :);
            mean_trace = squeeze(mean(diff_traces, 4, 'omitnan')); 
            se_trace = squeeze(std(diff_traces, [], 4, 'omitnan')./sqrt(size(saved_traces_all, 3)));
        else
            mean_trace = squeeze(mean(saved_traces_all(1:(end - 5), 2, tr_type, :), 4, 'omitnan')); 
            se_trace = squeeze(std(saved_traces_all(1:(end - 5), 2, tr_type, :), [], 4, 'omitnan')./sqrt(size(saved_traces_all, 3)));
        end
        
        plt_h = shadedErrorBar([], mean_trace, se_trace, {'Color', color_vecs(tr_type, :)}, 1);
        set(plt_h.edge(:), 'Color', 'none')
        
        write_data_cols = [write_data_cols, mean_trace, se_trace];  %logging data to write to file
        
        if plot_diff_traces == 1
            break
        else
            hold on
        end
    end
    hold off
    ylabel('unpaired odor responses (dF/F)')
    set_xlabels_time(2, frame_time, 10)
    ax_vals = axis;
    ax_vals(2) = 290;
    ax_vals(4) = 3;
    ax_vals(3) = -0.5;
    axis(ax_vals);
    fig_wrapup(2, [], [75, 90], 0.6)
    if set_list_type == 0
        add_stim_bar(2, stim_frs_bar, unpaired_color);
    elseif set_list_type == 1
        add_stim_bar(2, stim_frs_bar, [paired_color; unpaired_color]);
    elseif set_list_type == 2
        add_stim_bar(2, stim_frs_bar, [EL_color; unpaired_color]);
    elseif set_list_type == 3
        add_stim_bar(2, stim_frs_bar, [unpaired_color; EL_color]);
    elseif set_list_type == 4   %generalization datsets
        add_stim_bar(2, stim_frs_bar, [unpaired_color; EL_color]);  
    end
    
    
    if set_list_type == 0
        header_row = [{'Ap_pre_mean'}, {'Ap_pre_se'}, {'Ap_post_mean'}, {'Ap_post_se'}];
        xls_path = [paper_save_dir_simple,  'Ap_traces_', set_type_name, '.xls'];
    elseif set_list_type == 1
        header_row = [{'AAp_pre_mean'}, {'AAp_pre_se'}, {'AAp_post_mean'}, {'AAp_post_se'}];
        xls_path = [paper_save_dir,  'AAp_traces_', set_type_name, '.xls'];
    elseif set_list_type == 2
        header_row = [{'BAp_pre_mean'}, {'BAp_pre_se'}, {'BAp_post_mean'}, {'BAp_post_se'}];
        xls_path = [paper_save_dir,  'BAp_traces_', set_type_name, '.xls'];
    elseif set_list_type == 3
        header_row = [{'ApB_pre_mean'}, {'ApB_pre_se'}, {'ApB_post_mean'}, {'ApB_post_se'}];
        xls_path = [paper_save_dir,  'ApB_traces_', set_type_name, '.xls'];
    else
    end
    %writing data behind plot to file
    [c] = write_xls_header(header_row, write_data_cols, xls_path);
    write_data_cols = [];
     
    
    %EL
    figure(3)
    for tr_type = 1:2
        mean_trace = squeeze(mean(saved_traces_all(1:(end - 20), 3, tr_type, :), 4, 'omitnan')); 
        se_trace = squeeze(std(saved_traces_all(1:(end - 20), 3, tr_type, :), [], 4, 'omitnan')./sqrt(size(saved_traces_all, 3)));
        col_mult = 0.6.^(tr_type - 1);
        plt_h = shadedErrorBar([], mean_trace, se_trace, {'Color', color_vecs(tr_type, :)}, 1);
        set(plt_h.edge(:), 'Color', 'none')
        hold on
        
         write_data_cols = [write_data_cols, mean_trace, se_trace];  %logging data to write to file
    end
    hold off
    ylabel('EL responses (dF/F)')
    set_xlabels_time(3, frame_time, 10)
    ax_vals = axis;
    ax_vals(2) = 290;
    ax_vals(4) = 6;
    ax_vals(3) = -0.5;
    axis(ax_vals);
    fig_wrapup(3, [], [75, 90], 0.6)
    if set_list_type ~= 4
        add_stim_bar(3, stim_frs_bar_EL(1, :), EL_color);
    elseif set_list_type == 4
        add_stim_bar(3, stim_frs_bar_EL(1, :), paired_color);
    else
    end
    
     if set_list_type < 2
        header_row = [{'B_pre_mean'}, {'B_pre_se'}, {'B_post_mean'}, {'B_post_se'}];
        xls_path = [paper_save_dir,  'B_traces_', set_type_name, '.xls'];
     else
     end
    %writing data behind plot to file
    [c] = write_xls_header(header_row, write_data_cols, xls_path);
    write_data_cols = [];
    
    
     
    %quantification and statistical testing
    
    %comparing pre with post
    %reshaping resp size matrix for plotting function
    resp_mat_small = [];
    if set_list_type <= 1
        n_ods = 3;
    elseif set_list_type > 1
        n_ods = 2;
    else
    end
    
    for odor_n = 1:n_ods
        for resp_type = 1:2
            curr_resps = squeeze(resp_mat_all(resp_type, odor_n, :));
            resp_mat_small = [resp_mat_small, curr_resps];
        end
    end
    
    paired_multiplier = 1;
    if set_list_type <=2            %including EL first
        marker_colors = [paired_color; paired_color.*paired_multiplier; unpaired_color; unpaired_color.*paired_multiplier; EL_color; EL_color.*paired_multiplier];
        marker_colors = marker_colors(1:(n_ods.*2), :);
    elseif set_list_type == 3       %EL second
        marker_colors = repmat(EL_color, 6, 1);       %using grey markers because EL is always the second pulse
        marker_colors = marker_colors(1:(n_ods.*2), :);
    elseif set_list_type == 4       %generalization experiment
        marker_colors = [unpaired_color; unpaired_color.*paired_multiplier; EL_color; EL_color.*paired_multiplier];       
        marker_colors = marker_colors(1:(n_ods.*2), :);
    else
    end
        
    line_colors = zeros(size(marker_colors, 1), size(marker_colors, 2)) + 0.7;
    col_pairs = [1, 2; 3, 4; 5, 6];
    col_pairs = col_pairs(1:n_ods, :);
    xlabels = [{'prd pre'}, {'prd post'}, {'unprd pre'}, {'unprd post'}, {'EL pre'}, {'EL post'}];
    xlabels = xlabels(1:(n_ods.*2));
    figure(4)
    fig_h = scattered_dot_plot_ttest(resp_mat_small, 4, 0.6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 1, 0.05, 0, 1, 'force_mean');
    ylabel('response size (dF/F)');

    fig_wrapup(fig_h, [], [100, 120], 0.6);
    
    write_data_cols = resp_mat_small;
    if set_list_type == 0
        header_row = [{'A_pre'}, {'A_post'}, {'Ap_pre'}, {'Ap_post'}, {'B_pre'}, {'B_post'}];
        xls_path = [paper_save_dir_simple,  'resp_sizes_', set_type_name, '.xls'];
    elseif set_list_type == 1
        header_row = [{'ApA_pre'}, {'ApA_post'}, {'AAp_pre'}, {'AAp_post'}, {'B_pre'}, {'B_post'}];
        xls_path = [paper_save_dir,  'resp_sizes_', set_type_name, '.xls'];
    elseif set_list_type == 2
        header_row = [{'BA_pre'}, {'BA_post'}, {'BAp_pre'}, {'BAp_post'}];
        xls_path = [paper_save_dir,  'resp_sizes_', set_type_name, '.xls'];
    elseif set_list_type == 3
        header_row = [{'AB_pre'}, {'AB_post'}, {'ApB_pre'}, {'ApB_post'}];
        xls_path = [paper_save_dir,  'resp_sizes_', set_type_name, '.xls'];
    else
    end
    
    %writing data behind plot to file
    [c] = write_xls_header(header_row, write_data_cols, xls_path);
    write_data_cols = [];
      
    %computing post means and ses as a percentage of pre means
    mean_vals = mean(resp_mat_small, 1, 'omitnan');
    SE_vals = std(resp_mat_small, [], 1, 'omitnan')./sqrt(size(resp_mat_small, 1));
    mean_percent_vec = [];
    se_percent_vec = [];
    for col_ni = 1:(size(resp_mat_small, 2)./2)
        col_n = col_ni.*2;  %only computing for post pairing columns
        curr_premean = mean_vals((col_n - 1));
        mean_percent_vec = [mean_percent_vec; (mean_vals(col_n)./curr_premean)];
        se_percent_vec = [se_percent_vec; (SE_vals(col_n)./curr_premean)];
    end
    display(mean_percent_vec)
    display(se_percent_vec)
    
    
    %comparing paired odor resps with unpaired odor resps
    resp_mat_small_sw = [resp_mat_small(:, 1), resp_mat_small(:, 3), resp_mat_small(:, 2), resp_mat_small(:, 4)];
    marker_colors_sw = [marker_colors(1, :); marker_colors(3, :); marker_colors(2, :); marker_colors(4, :)];
    col_pairs = [1, 2; 3, 4];
    xlabels = [{'prd pre'}, {'unprd pre'}, {'prd post'}, {'unprd post'}];
    figure(5)
    fig_h = scattered_dot_plot_ttest(resp_mat_small_sw, 5, .6, 1, 4, marker_colors_sw, 1, col_pairs, line_colors(1:4, :), xlabels, 2, mean_color, 1, 0.05, 0, 1, 'force_mean');
    ylabel('response size (dF/F)');
    fig_wrapup(fig_h, [], [100, 80], 0.6);
    
     
    
    %Plotting a contrast score (pulse2 - pulse1 response) for transition
    %trials
    contrast_scores = ctrst_mat_all_pulse2 - resp_mat_all_pulse1;
    cscore_mat = [squeeze(contrast_scores(1, 1, :)), squeeze(contrast_scores(1, 2, :)), squeeze(contrast_scores(2, 1, :)), squeeze(contrast_scores(2, 2, :))];
    paired_multiplier = 1;
    if set_list_type ~= 4
        marker_colors = [ paired_color; unpaired_color; paired_color.*paired_multiplier;  unpaired_color.*paired_multiplier];
    elseif set_list_type == 4   %for generalization experiment
        marker_colors = [ unpaired_color; EL_color; unpaired_color.*paired_multiplier;  EL_color.*paired_multiplier];
    end
    line_colors = repmat([0.6, 0.6, 0.6], 4, 1);
    col_pairs = [1, 2; 3, 4];
    xlabels = [{'prd pre'}, {'unprd pre'}, {'prd post'}, {'unprd post'}];
    figure(7)
    fig_h = scattered_dot_plot_ttest(cscore_mat, 7, 0.6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 1, 0.05, 0, 1, 'force_mean');
    ylabel('contrast at transition (dF/F)');
    fig_wrapup(fig_h, [], [100, 80], 0.6);
    
    write_data_cols = cscore_mat;
    if set_list_type == 0
        header_row = [{'A_pre'}, {'A_post'}, {'Ap_pre'}, {'Ap_post'}];
        xls_path = [paper_save_dir_sfig,  'contrast_scores_', set_type_name, '.xls'];
    elseif set_list_type == 1
        header_row = [{'ApA_pre'}, {'ApA_post'}, {'AAp_pre'}, {'AAp_post'}];
        xls_path = [paper_save_dir_sfig,  'contrast_scores_', set_type_name, '.xls'];
    elseif set_list_type == 2
        header_row = [{'BA_pre'}, {'BA_post'}, {'BAp_pre'}, {'BAp_post'}];
        xls_path = [paper_save_dir_sfig,  'contrast_scores_', set_type_name, '.xls'];
    elseif set_list_type == 3
        header_row = [{'AB_pre'}, {'AB_post'}, {'ApB_pre'}, {'ApB_post'}];
        xls_path = [paper_save_dir_sfig,  'contrast_scores_', set_type_name, '.xls'];
    else
    end
    
    %writing data behind plot to file
    [c] = write_xls_header(header_row, write_data_cols, xls_path);
    write_data_cols = [];
    
        
    %plotting individual traces
    od_type_names = [{'paired '}, {'unpaired '}];
    stim_type_names = [{'pre'}, {'post'}];
    color_vecs = [paired_color; unpaired_color];

    %sorting single traces by total activity in the unpaired, post response
    %dataset
    curr_mat = squeeze(saved_traces_all(:, 2, 2, :))';
    tot_act = sum(curr_mat, 2, 'omitnan');
    
    for od_type = 1:2   %paired-unpaired loop
        for stim_type = 1:2 %pre-post loop
            curr_mat = squeeze(saved_traces_all(:, od_type, stim_type, :))';
            curr_PID_mat = squeeze(saved_PID_traces_all(:, od_type, stim_type, :))';
                       
            fig_h = figure('Name', [od_type_names{od_type}, stim_type_names{stim_type}, ' single traces']);
            curr_mat = [tot_act, curr_mat];
            curr_mat = sortrows(curr_mat);
            curr_mat(:, 1) = [];
            if size(curr_mat, 2) > 490
                curr_mat(:, 490:end) = [];
            else
            end
            cascade_plot(fig_h, curr_mat', [0.6, 0.6, 0.6], 1, 0, .3, 3, 3, 0, 1);
            set_xlabels_time(fig_h, frame_time, 10);
            other_type = [1, 2];
            other_type(other_type == od_type) = [];
            %fig_wrapup(fig_h, []);
            
            if set_list_type == 1
                add_stim_bar(fig_h, stim_frs_bar, [color_vecs(other_type, :); color_vecs(od_type, :)] );
            elseif set_list_type == 0
                add_stim_bar(fig_h, stim_frs_bar, [color_vecs(od_type, :)] );
            elseif set_list_type == 2
                add_stim_bar(fig_h, stim_frs_bar, [EL_color; color_vecs(od_type, :)]);
            elseif set_list_type == 3
                add_stim_bar(fig_h, stim_frs_bar, [color_vecs(od_type, :); EL_color]);
            elseif set_list_type == 4   %generalization expts
                if od_type == 1
                    add_stim_bar(fig_h, stim_frs_bar, [EL_color; color_vecs(2, :)]);
                elseif od_type == 2
                    add_stim_bar(fig_h, stim_frs_bar, [color_vecs(2, :); EL_color]);
                else
                end
            end
            
        end
    end
    
    
end