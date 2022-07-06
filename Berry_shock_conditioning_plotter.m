clear all
close all


%MBONA3 imaging
dataset_list_paths = [...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_WTalpha3_fedshock_PAplus.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_WTalpha3_fedshock_BAplus.xls'};...
                      {'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_WTalpha3_fedshock_unpairedctrl.xls'};...
                      ];
                  

%MBONG2A'1 imaging (1x pairing)                 
%dataset_list_paths = [...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_WTalpha3_fedshock_PAplus.xls'};...       %NOTE: this is a DUMMY set list to be replaced with real list once acquired
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_WTG2Ap1_fedshock_BAplus_1x.xls'};...
                      %{'C:\Data\Code\general_code_old\data_folder_lists\Janelia\Berry_handover_WTG2Ap1_fedshock_unpairedctrl_1x.xls'};...
    
 %                     ];
            
suppress_plots = 1;
plotting_quant_no_filt = 0;     %1 - only unfiltered traces used for all analysis and plotting - traces included. 0 - filtered traces used for everything.

% [del, odor_names1] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
% [del, odor_names2] = xlsread('C:\Data\Code\general_code_old\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);
% odor_names2{3} = 'Butyl acetate';

cell_n = 1;

global color_vec;                
a = colormap('bone');
global greymap
greymap = flipud(a);
fly_n = 0;
script_name = mfilename;

paired_color = [0.851, 0.3725, 0.0078];
unpaired_color = [0.1059, 0.6196, 0.4667];
EL_color = [0.4588, 0.4392, 0.7020];
PA_color = [0.2667, 0.9569, 0.9255].*0.8;
BA_color = [0.5549, 0.9686, 0.433].*0.8;
mean_color = ([0, 49, 152]./256).*1.5;

paired_color = [0,136,55]./256;
unpaired_color = [166,219,160]./256;
EL_color = [123,50,148]./256;
mean_color = [0, 0, 0];
%mean_color = [0.8, 0.4, 0.4];

y_ax_traces = 0.8;
y_ax_fit_traces = 0.6;

an_save_path = 'C:\Data\Data\Analysed_data\Analysis_results\PaBaEl_Gamma2\';
paper_save_dir = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\Fig4_MBON_transitions_fig\';
paper_save_dir_sfig = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\SFig3_4_MBON_transitions_fig\';
force_resave = 1;

n_sec = 2;      %width of moving integration window in s

integ_win_s = 5;        %width of pulse integration window for response quantification
saved_traces_all = [];
list_eg = dataset_list_paths{1};
list_set_namei = findstr(list_eg, 'WT');
list_set_name = list_eg((list_set_namei + 2):(list_set_namei + 7));

if exist(['C:\Data\Data\Analysed_data\Analysis_results\fine_discr_shock_cond\', list_set_name,'\saved_traces_all.mat']) ~= 2
    for list_n = 1:size(dataset_list_paths, 1)

        curr_dir_list_path = dataset_list_paths{list_n, 1};
        [del, dir_list] = xlsread(curr_dir_list_path, 1);        %list of Suite2P results directories
        n_dirs = size(dir_list, 1);
        dataset_list_name = findstr(curr_dir_list_path, 'list_');
        dataset_list_name = curr_dir_list_path((dataset_list_name + 5):(end - 4));

        if isempty(findstr(curr_dir_list_path, 'BA')) == 0
            paired_od_n_olf2 = 3;
            unpaired_od_n_olf2 = 1;
            set_list_type = 1;      %paired list
        elseif isempty(findstr(curr_dir_list_path, 'PA')) == 0
            paired_od_n_olf2 = 1;
            unpaired_od_n_olf2 = 3;
            set_list_type = 1;      %paired list
        elseif  isempty(findstr(curr_dir_list_path, 'unpaired')) == 0
            set_list_type = 0;  %unpaired list
            paired_od_n_olf2 = 1;   %this is just a dummy, to determine indexing
            unpaired_od_n_olf2 = 3;
        end

        saved_traces_list = [];
        saved_PID_traces_all = [];
        %loop to go through all experiment datasets listed in list file
        for dir_n = 1:n_dirs
            fly_n = fly_n + 1;

            curr_dir = [dir_list{dir_n, 1}, '\'];
            curr_dir = manage_base_paths(curr_dir, 2)
            saved_traces = [];


            tif_times = load([curr_dir, 'tif_time_stamps.mat']);           %reading in time stamps for each tif file recorded by raw_data_extracter
            tif_times = tif_times.time_stamps;

            cd(curr_dir);
            tif_name = dir('*.tif');
            if isempty(tif_name) == 1
                keyboard
                continue
            else
            end
            stack_obj = ScanImageTiffReader([curr_dir, tif_name(1).name]);
            [frame_time, zoom, n_chans, PMT_offsets] = SI_tif_info(stack_obj);

            [stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig, PID_traces] = load_params_trains_modular(curr_dir, tif_times, frame_time);    %reading in trial stimulus parameters after matching time stamps to F traces

            %looking up odor numbers across olfactometers based on current dataset's odor name list
            odor_names1 = stim_mat.odourNames;
            odor_names2 = stim_mat.odourNames_olf2;
            paired_od_n_olf1 = od_name_lookup(odor_names1, odor_names2{paired_od_n_olf2});
            unpaired_od_n_olf1 = od_name_lookup(odor_names1, odor_names2{unpaired_od_n_olf2});

            paired_color = color_vec(2, :);
            unpaired_color = color_vec(1, :);
            EL_color = color_vec(3, :);
            %mean_color = [0.5, 0.83, 0.98];
            %mean_color = [149, 200, 216]./256;
            mean_color = [0, 0, 0];
            %mean_color = ([0, 49, 152]./256).*1.5;
           
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

            %if weird EL trial1 exists, getting rid of it
            if stim_mat_simple(1, od_olf1_col_n) == 11
                stim_mat(1) = [];
                stim_mat_simple(1, :) = [];
                raw_data_mat(:, :, 1) = [];
            else
            end

            good_tr_list = [1, 2, 3, 6, 7, 8];

            n_cells = size(raw_data_mat, 2);

            %calculating dF/F traces from raw data
            filt_time = 0.2;            %in s, the time window for boxcar filter for generating filtered traces
            [dff_data_mat, dff_data_mat_f] = cal_dff_traces_res(raw_data_mat, stim_mat, frame_time, filt_time, curr_dir, good_tr_list);
            if plotting_quant_no_filt == 1
                dff_data_mat_f = dff_data_mat;
            else
            end

    %         if size(dff_data_mat, 2) > 1
    %             dff_data_mat = mean(dff_data_mat, 2, 'omitnan');
    %             dff_data_mat_f = mean(dff_data_mat_f, 2, 'omitnan');
    %         else
    %         end
            dff_data_mat = dff_data_mat(:, cell_n, :);
            dff_data_mat_f  = dff_data_mat_f(:, cell_n, :);
            del = find(dff_data_mat_f < -1);
            dff_data_mat_f(del) = -1;       %forcing crazy values to sane ones

            y_ax_lim = [];
            plot_means = 1;

            %plotting, quantifying resps

            %1. paired odor response
            curr_tr = find(stim_mat_simple(1:3, od_col_ns(2) ) == paired_od_n_olf2);
            curr_trace = dff_data_mat_f(:, 1, curr_tr);
            saved_traces = pad_n_concatenate(saved_traces, curr_trace, 2, nan);

            %2. paired odor response
            curr_tr = find(stim_mat_simple(1:3, od_col_ns(2) ) == unpaired_od_n_olf2);
            curr_trace = dff_data_mat_f(:, 1, curr_tr);
            saved_traces = pad_n_concatenate(saved_traces, curr_trace, 2, nan);

            %3. EL odor response
            curr_tr = find(stim_mat_simple(1:3, od_col_ns(2) ) == 4);
            curr_trace = dff_data_mat_f(:, 1, curr_tr);
            saved_traces = pad_n_concatenate(saved_traces, curr_trace, 2, nan);


            %logging data across flies
            saved_traces_list = pad_n_concatenate(saved_traces_list, saved_traces, 3, nan);

            %plotting traces
            stim_frs = compute_stim_frs_modular(stim_mat, 1, frame_time);
            if suppress_plots == 0                
                stim_frs_bar = [stim_frs{1}; stim_frs{2}];
                %paired_od
                figure(1)
                plot(saved_traces(:, 1), 'Color', paired_color, 'lineWidth', 2);
                ylabel('paired odor responses (dF/F)')
                set_xlabels_time(1, frame_time, 10)
                ax_vals = axis;
                ax_vals(4) = 6;
                ax_vals(3) = -2;
                axis(ax_vals);
                fig_wrapup(1, [], [50, 60])
                add_stim_bar(1, stim_frs_bar, [unpaired_color; paired_color]);


                %unpaired_od
                figure(2)
                plot(saved_traces(:, 2), 'Color', unpaired_color, 'lineWidth', 2);
                ylabel('unpaired odor responses (dF/F)')
                set_xlabels_time(2, frame_time, 10)
                ax_vals = axis;
                ax_vals(4) = 6;
                ax_vals(3) = -2;
                axis(ax_vals);
                fig_wrapup(2, [], [50, 60])
                add_stim_bar(2, stim_frs_bar, [paired_color; unpaired_color]);


                %EL
                figure(3)
                plot(saved_traces(:, 3), 'Color', EL_color, 'lineWidth', 2);
                ylabel('EL responses (dF/F)')
                set_xlabels_time(3, frame_time, 10)
                ax_vals = axis;
                ax_vals(4) = 6;
                ax_vals(3) = -2;
                axis(ax_vals);
                fig_wrapup(3, [], [50, 60]) 
                add_stim_bar(3, stim_frs{1}, EL_color);

                 keyboard
                 close figure 1
                 close figure 2
                 close figure 3
            else
            end

        end
        saved_traces_all = pad_n_concatenate(saved_traces_all, saved_traces_list, 4, nan);
    end
    mkdir(['C:\Data\Data\Analysed_data\Analysis_results\fine_discr_shock_cond\', list_set_name]);
    save(['C:\Data\Data\Analysed_data\Analysis_results\fine_discr_shock_cond\', list_set_name,'\saved_traces_all.mat'], 'saved_traces_all');
    save(['C:\Data\Data\Analysed_data\Analysis_results\fine_discr_shock_cond\', list_set_name,'\stim_frs.mat'], 'stim_frs');
    save(['C:\Data\Data\Analysed_data\Analysis_results\fine_discr_shock_cond\', list_set_name,'\frame_time.mat'], 'frame_time');
else
    saved_traces_all = load(['C:\Data\Data\Analysed_data\Analysis_results\fine_discr_shock_cond\', list_set_name,'\saved_traces_all.mat']);
    saved_traces_all = saved_traces_all.saved_traces_all;
    
    stim_frs = load(['C:\Data\Data\Analysed_data\Analysis_results\fine_discr_shock_cond\', list_set_name,'\stim_frs.mat']);
    stim_frs = stim_frs.stim_frs;
    
    frame_time = load(['C:\Data\Data\Analysed_data\Analysis_results\fine_discr_shock_cond\', list_set_name,'\frame_time.mat']);
    frame_time = frame_time.frame_time;
end

%saved_traces_all dim1: time, dim2: [paired_od, unpaired_od, EL], dim3: fly_n, dim4: [BApaired, PApaired, unpairedctrl]

%computing number of real flies for each dataset list type
n_flies = size(saved_traces_all, 3) - squeeze(sum(isnan(saved_traces_all(1, 1, :, :)), 3));

%plotting mean response traces
%1. Plotting by pairing status

%computing mean and se traces for paired datasets
mean_paired_trace = squeeze(mean(mean(saved_traces_all(:, 1, :, 1:2), 3, 'omitnan'), 4, 'omitnan'));
mean_unpaired_trace = squeeze(mean(mean(saved_traces_all(:, 2, :, 1:2), 3, 'omitnan'), 4, 'omitnan'));
mean_EL_trace = squeeze(mean(mean(saved_traces_all(:, 3, :, 1:2), 3, 'omitnan'), 4, 'omitnan'));

var_paired_trace = squeeze(mean(var(saved_traces_all(:, 1, :, 1:2), 0, 3, 'omitnan'), 4, 'omitnan'));
se_paired_trace = sqrt(var_paired_trace)./sqrt(n_flies(1) + n_flies(2));
var_unpaired_trace = squeeze(mean(var(saved_traces_all(:, 2, :, 1:2), 0, 3, 'omitnan'), 4, 'omitnan'));
se_unpaired_trace = sqrt(var_unpaired_trace)./sqrt(n_flies(1) + n_flies(2));
var_EL_trace = squeeze(mean(var(saved_traces_all(:, 3, :, 1:2), 0, 3, 'omitnan'), 4, 'omitnan'));
se_EL_trace = sqrt(var_EL_trace)./sqrt(n_flies(1) + n_flies(2));

%computing mean and se traces for unpaired datasets
mean_sim_od_trace_ct = squeeze(mean(mean(saved_traces_all(:, 1:2, :, 3), 3, 'omitnan'), 2, 'omitnan'));
mean_EL_trace_ct = squeeze(mean(saved_traces_all(:, 3, :, 3), 3, 'omitnan'));

var_sim_od_trace_ct = squeeze(mean(var(saved_traces_all(:, 1:2, :, 3), 0, 3, 'omitnan'), 2, 'omitnan'));
se_sim_od_trace_ct = sqrt(var_sim_od_trace_ct)./sqrt(n_flies(1) + n_flies(2));
var_EL_trace_ct = squeeze(var(saved_traces_all(:, 3, :, 3), 0, 3, 'omitnan'));
se_EL_trace_ct = sqrt(var_EL_trace_ct)./sqrt(n_flies(3));


%plotting
stim_frs_bar = [stim_frs{1}; stim_frs{2}];
figure(1)
plt_h = shadedErrorBar([], mean_sim_od_trace_ct, se_sim_od_trace_ct, {'Color', [0.6, 0.6, 0.6]}, 1);
set(plt_h.edge(:), 'Color', 'none')
hold on
plt_h = shadedErrorBar([], mean_paired_trace, se_paired_trace, {'Color', [0, 0, 0]}, 1);
set(plt_h.edge(:), 'Color', 'none')
hold off
ylabel('paired odor responses (dF/F)')
set_xlabels_time(1, frame_time, 10)
ax_vals = axis;
ax_vals(2) = 290;
ax_vals(4) = 3.3;
ax_vals(3) = -0.5;
axis(ax_vals);
fig_wrapup(1, [], [75, 90], 0.6)
add_stim_bar(1, stim_frs_bar, [unpaired_color; paired_color]);

header_row = [{'ApA_pre_mean'}, {'ApA_pre_se'}, {'ApA_post_mean'}, {'ApA_post_se'}];
write_data_cols = [];
write_data_cols = mean_sim_od_trace_ct;
write_data_cols = pad_n_concatenate(write_data_cols, se_sim_od_trace_ct, 2, nan);
write_data_cols = pad_n_concatenate(write_data_cols, mean_paired_trace, 2, nan);
write_data_cols = pad_n_concatenate(write_data_cols, se_paired_trace, 2, nan);

%writing data behind plot to file
xls_path = [paper_save_dir,  'ApA_traces_transition_MBONA3.xls'];
[c] = write_xls_header(header_row, write_data_cols, xls_path);
write_data_cols = [];




figure(2)
plt_h = shadedErrorBar([], mean_sim_od_trace_ct, se_sim_od_trace_ct, {'Color', [0.6, 0.6, 0.6]}, 1);
set(plt_h.edge(:), 'Color', 'none')
hold on
plt_h = shadedErrorBar([], mean_unpaired_trace, se_unpaired_trace, {'Color', [0, 0, 0]}, 1);
set(plt_h.edge(:), 'Color', 'none')
hold off
ylabel('unpaired odor responses (dF/F)')
set_xlabels_time(2, frame_time, 10)
ax_vals = axis;
ax_vals(2) = 290;
ax_vals(4) = 3.3;
ax_vals(3) = -0.5;
axis(ax_vals);
fig_wrapup(2, [], [75, 90], 0.6)
add_stim_bar(2, stim_frs_bar, [paired_color; unpaired_color]);

header_row = [{'AAp_pre_mean'}, {'AAp_pre_se'}, {'AAp_post_mean'}, {'AAp_post_se'}];
write_data_cols = [];
write_data_cols = pad_n_concatenate(write_data_cols, mean_sim_od_trace_ct, 2, nan);
write_data_cols = pad_n_concatenate(write_data_cols, se_sim_od_trace_ct, 2, nan);
write_data_cols = pad_n_concatenate(write_data_cols, mean_unpaired_trace, 2, nan);
write_data_cols = pad_n_concatenate(write_data_cols, se_unpaired_trace, 2, nan);

%writing data behind plot to file
xls_path = [paper_save_dir,  'AAp_traces_transition_MBONA3.xls'];
[c] = write_xls_header(header_row, write_data_cols, xls_path);
write_data_cols = [];


figure(3)
plt_h = shadedErrorBar([], mean_EL_trace_ct, se_EL_trace_ct, {'Color', [0.6, 0.6, 0.6]}, 1);
set(plt_h.edge(:), 'Color', 'none')
hold on
plt_h = shadedErrorBar([], mean_EL_trace, se_EL_trace, {'Color', [0, 0, 0]}, 1);
set(plt_h.edge(:), 'Color', 'none')
hold off
ylabel('unpaired odor responses (dF/F)')
set_xlabels_time(3, frame_time, 10)
ax_vals = axis;
ax_vals(2) = 290;
ax_vals(4) = 4;
ax_vals(3) = -0.5;
axis(ax_vals);
fig_wrapup(3, [], [75, 90], 0.6)
add_stim_bar(3, stim_frs{1}, [EL_color]);
header_row = [{'B_pre_mean'}, {'B_pre_se'}, {'B_post_mean'}, {'B_post_se'}];
write_data_cols = [];
write_data_cols = pad_n_concatenate(write_data_cols, mean_EL_trace_ct, 2, nan);
write_data_cols = pad_n_concatenate(write_data_cols, se_EL_trace_ct, 2, nan);
write_data_cols = pad_n_concatenate(write_data_cols, mean_EL_trace, 2, nan);
write_data_cols = pad_n_concatenate(write_data_cols, se_EL_trace, 2, nan);

%writing data behind plot to file
xls_path = [paper_save_dir,  'B_traces_transitions_MBONA3.xls'];
[c] = write_xls_header(header_row, write_data_cols, xls_path);
write_data_cols = [];


%2. Plotting by odor identity
%PA
mean_PA_traces1 = mean(saved_traces_all(:, 1, :, 1), 3, 'omitnan');     %PA paired
mean_PA_traces2 = mean(saved_traces_all(:, 2, :, 2), 3, 'omitnan');     %PA not paired
%mean_PA_traces = mean([mean_PA_traces1, mean_PA_traces2], 2, 'omitnan');
var_PA_traces1 = var(saved_traces_all(:, 1, :, 1), 0, 3, 'omitnan');
se_PA_traces1 = sqrt(var_PA_traces1)./sqrt(n_flies(1));
var_PA_traces2 = var(saved_traces_all(:, 1, :, 1), 0, 3, 'omitnan');
se_PA_traces2 = sqrt(var_PA_traces2)./sqrt(n_flies(2));
% var_PA_traces = mean([var_PA_traces1, var_PA_traces2], 2, 'omitnan');
% se_PA_traces = sqrt(var_PA_traces)./sqrt(n_flies(1) + n_flies(2));

mean_PA_traces_ct = mean(saved_traces_all(:, 1, :, 3), 3, 'omitnan');   %note, for indexing the ctrl datasets, PA was treated as the 'paired odor'
var_PA_traces_ct = var(saved_traces_all(:, 1, :, 3), 0, 3, 'omitnan');
se_PA_traces_ct = sqrt(var_PA_traces_ct)./sqrt(n_flies(3));


%BA
mean_BA_traces1 = mean(saved_traces_all(:, 1, :, 2), 3, 'omitnan');     %BA paired
mean_BA_traces2 = mean(saved_traces_all(:, 2, :, 1), 3, 'omitnan');     %BA not paired
%mean_BA_traces = mean([mean_BA_traces1, mean_BA_traces2], 2, 'omitnan');
var_BA_traces1 = var(saved_traces_all(:, 1, :, 2), 0, 3, 'omitnan');
se_BA_traces1 = sqrt(var_BA_traces1)./sqrt(n_flies(2));
var_BA_traces2 = var(saved_traces_all(:, 2, :, 1), 0, 3, 'omitnan');
se_BA_traces2 = sqrt(var_BA_traces2)./sqrt(n_flies(1));
% var_BA_traces = mean([var_BA_traces1, var_BA_traces2], 2, 'omitnan');
% se_BA_traces = sqrt(var_BA_traces)./sqrt(n_flies(1) + n_flies(2));

mean_BA_traces_ct = mean(saved_traces_all(:, 2, :, 3), 3, 'omitnan');   %note, for indexing the ctrl datasets, PA was treated as the 'paired odor'
var_BA_traces_ct = var(saved_traces_all(:, 2, :, 3), 0, 3, 'omitnan');
se_BA_traces_ct = sqrt(var_BA_traces_ct)./sqrt(n_flies(3));


%EL
mean_EL_traces1 = mean(saved_traces_all(:, 3, :, 1), 3, 'omitnan');
mean_EL_traces2 = mean(saved_traces_all(:, 3, :, 2), 3, 'omitnan');
mean_EL_traces = mean([mean_EL_traces1, mean_EL_traces2], 2, 'omitnan');
var_EL_traces1 = var(saved_traces_all(:, 3, :, 1), 0, 3, 'omitnan');
var_EL_traces2 = var(saved_traces_all(:, 3, :, 1), 0, 3, 'omitnan');
var_EL_traces = mean([var_EL_traces1, var_EL_traces2], 2, 'omitnan');
se_EL_traces = sqrt(var_EL_traces)./sqrt(n_flies(1) + n_flies(2));

mean_EL_traces_ct = mean(saved_traces_all(:, 3, :, 3), 3, 'omitnan');   %note, for indexing the ctrl datasets, PA was treated as the 'paired odor'
var_EL_traces_ct = var(saved_traces_all(:, 3, :, 3), 0, 3, 'omitnan');
se_EL_traces_ct = sqrt(var_EL_traces_ct)./sqrt(n_flies(3));


%Plotting
%PA
figure(4)
plt_h = shadedErrorBar([], mean_PA_traces_ct, se_PA_traces_ct, {'Color', [0.6, 0.6, 0.6]}, 1);
set(plt_h.edge(:), 'Color', 'none')
hold on
plt_h = shadedErrorBar([], mean_PA_traces1, se_PA_traces1, {'Color', [0, 0, 0]}, 1);
set(plt_h.edge(:), 'Color', 'none')
plt_h = shadedErrorBar([], mean_PA_traces2, se_PA_traces2, {'--', 'Color', [0, 0, 0]}, 1);
set(plt_h.edge(:), 'Color', 'none')
hold off
ylabel('PA responses (dF/F)')
set_xlabels_time(4, frame_time, 10)
ax_vals = axis;
ax_vals(2) = 290;
ax_vals(4) = 4;
ax_vals(3) = -0.5;
axis(ax_vals);
fig_wrapup(4, [], [75, 90], 0.6)
add_stim_bar(4, stim_frs_bar, [BA_color; PA_color]);

%BA
figure(5)
plt_h = shadedErrorBar([], mean_BA_traces_ct, se_BA_traces_ct, {'Color', [0.6, 0.6, 0.6]}, 1);
set(plt_h.edge(:), 'Color', 'none')
hold on
plt_h = shadedErrorBar([], mean_BA_traces1, se_BA_traces1, {'Color', [0, 0, 0]}, 1);
set(plt_h.edge(:), 'Color', 'none')
plt_h = shadedErrorBar([], mean_BA_traces2, se_BA_traces2, {'--', 'Color', [0, 0, 0]}, 1);
set(plt_h.edge(:), 'Color', 'none')
hold off
ylabel('BA responses (dF/F)')
set_xlabels_time(5, frame_time, 10)
ax_vals = axis;
ax_vals(2) = 290;
ax_vals(4) = 4;
ax_vals(3) = -0.5;
axis(ax_vals);
fig_wrapup(5, [], [75, 90], 0.6)
add_stim_bar(5, stim_frs_bar, [PA_color; BA_color]);


%STATISTICAL TESTING
%computing response sizes for statistical testing
%defining frame integration windows to compute response sizes
win1 = stim_frs{1};
win2 = stim_frs{2};
win2_orig = win2;
win2(2) = win2(2) + round(3./frame_time);

saved_resps_all1 = squeeze(mean(saved_traces_all(win1(1):win1(2), :, :, :), 1, 'omitnan'));      %mean responses during pulse1 time
saved_resps_all2 = squeeze(mean(saved_traces_all(win2(1):win2(2), :, :, :), 1, 'omitnan'));      %mean responses during pulse2 time

ctrst_sc_resps_all1 = squeeze(min(saved_traces_all(win1(1):win1(2), :, :, :), [], 1, 'omitnan'));      %mean responses during pulse1 time
ctrst_sc_resps_all2 = squeeze(max(saved_traces_all(win2_orig(1):win2_orig(2), :, :, :), [], 1, 'omitnan'));      %mean responses during pulse1 time

%saved_traces_all dim1: time, dim2: [paired_od, unpaired_od, EL], dim3: fly_n, dim4: [BApaired, PApaired, unpairedctrl]

%Plotting
%1. by pairing status
paired_resps = [saved_resps_all2(1, :, 1)'; saved_resps_all2(1, :, 2)'];
unpaired_resps = [saved_resps_all2(2, :, 1)'; saved_resps_all2(2, :, 2)'];

paired_ctrst_scores = [ctrst_sc_resps_all2(1, :, 1)'; ctrst_sc_resps_all2(1, :, 2)'] - [ctrst_sc_resps_all1(1, :, 1)'; ctrst_sc_resps_all1(1, :, 2)'];
unpaired_ctrst_scores = [ctrst_sc_resps_all2(2, :, 1)'; ctrst_sc_resps_all2(2, :, 2)'] - [ctrst_sc_resps_all1(2, :, 1)'; ctrst_sc_resps_all1(2, :, 2)'];

sim_resps_ctrl = [saved_resps_all2(1, :, 3)'; saved_resps_all2(2, :, 3)'];
EL_resps = [saved_resps_all1(3, :, 1)'; saved_resps_all1(3, :, 2)'];
EL_resps_ctrl = [saved_resps_all1(3, :, 3)'];

plot_mat = [sim_resps_ctrl, paired_resps, unpaired_resps];
plot_mat = pad_n_concatenate(plot_mat, EL_resps_ctrl, 2, nan);
plot_mat = pad_n_concatenate(plot_mat, EL_resps, 2, nan);

paired_multiplier = 1;
marker_colors = [[0.65, 0.65, 0.65]; paired_color.*paired_multiplier; unpaired_color.*paired_multiplier; [0.65, 0.65, 0.65]; EL_color.*paired_multiplier];
line_colors = [];
col_pairs = [];
xlabels = [{'ctrlsim'}, {'CSpl'}, {'CSmin'}, {'ctrldsim'}, {'dsim'}];
figure(7)
[fig_h, r_vecs_saved] = scattered_dot_plot_ttest(plot_mat, 7, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean');
%adding marker colors for each odor in col1 
hold on
curr_markersx = r_vecs_saved(1:n_flies(1), 1);
curr_markersy = plot_mat(1:n_flies(1), 1);              %Note: the ctrl column consists of ctrl PA resps with ctrl BA resps appended beneath
% plot(curr_markersx, curr_markersy, 'O', 'markerFaceColor', 'none', 'markerEdgeColor', [.65, .65, .65], 'markerSize', 4);
% curr_markersx = r_vecs_saved( (n_flies(1) + 1):end, 1);
% curr_markersy = plot_mat((n_flies(1) + 1):end, 1);      %Note: the ctrl column consists of ctrl PA resps with ctrl BA resps appended beneath
% plot(curr_markersx, curr_markersy, 'O', 'markerFaceColor', 'none', 'markerEdgeColor', [.65, .65, .65], 'markerSize', 4);
% mean_val = mean(plot_mat(:, 1), 'omitnan');
% se_val = std(plot_mat(:, 1), [], 'omitnan')./sqrt(size(plot_mat, 1));
% errorbar(1.75, mean_val, se_val, 'O', 'markerSize', 4, 'markerEdgeColor', mean_color, 'markerFaceColor', mean_color, 'Color', mean_color, 'lineWidth', 2)
hold off
ax_vals = axis;
ax_vals(3) = -0.5;
axis(ax_vals)
ylabel('response size (dF/F)');
fig_wrapup(fig_h, [], [100, 120], 0.6);
%testing
[ppaired, hpaired] = ranksum(plot_mat(:, 1), plot_mat(:, 2))
[punpaired, hunpaired] = ranksum(plot_mat(:, 1), plot_mat(:, 3))
[pEL, hEL] = ranksum(plot_mat(:, 4), plot_mat(:, 5))

p_corrected = bonf_holm([ppaired, punpaired, pEL], 0.05)


header_row = [{'ctrl_sim'}, {'A'}, {'Ap'}, {'ctrl_B'}, {'B'}];
write_data_cols = plot_mat;

%writing data behind plot to file
xls_path = [paper_save_dir,  'resp_sizes_transitions_MBONA3.xls'];
[c] = write_xls_header(header_row, write_data_cols, xls_path);
write_data_cols = [];



%comparing post A and A' with each other
figure(9)
[fig_h, r_vecs_saved] = scattered_dot_plot_ttest(plot_mat(:, 2:3), 9, .6, 1, 4, marker_colors(2:3, :), 1, col_pairs, line_colors, xlabels(2:3), 2, mean_color, 2, 0.05, 0, 1, 'force_mean');
ylabel('response size (dF/F)');
fig_wrapup(fig_h, [], [100, 60], 0.6);

%computing unpaired means and ses as a percentage of paired means
mean_vals = mean(plot_mat, 1, 'omitnan');
SE_vals = std(plot_mat, [], 1, 'omitnan')./sqrt(size(plot_mat, 1));
mean_percent_vec_post = [];
se_percent_vec_post = [];
pre_mean_vals = mean_vals([1, 4]);
for col_n = 1:(size(plot_mat, 2))
    if col_n < 4
        curr_premean = pre_mean_vals(1);
    else
        curr_premean = pre_mean_vals(2);
    end
      
    mean_percent_vec_post = [mean_percent_vec_post; (mean_vals(col_n)./curr_premean)];
    se_percent_vec_post = [se_percent_vec_post; (SE_vals(col_n)./curr_premean)];
end
display(mean_percent_vec_post)
display(se_percent_vec_post)



%plotting contrast scores
paired_multiplier = 1;
marker_colors = [paired_color; unpaired_color];
line_colors = [];
col_pairs = [];
xlabels = [{'CSpl'}, {'CSmin'}];
figure(8)
[fig_h, r_vecs_saved] = scattered_dot_plot_ttest([paired_ctrst_scores, unpaired_ctrst_scores], 8, 0.6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean');

ax_vals = axis;
ax_vals(3) = -0.5;
ax_vals(4) = 5;
axis(ax_vals)
ylabel('contrast score (dF/F)');
fig_wrapup(fig_h, [], [100, 60], 0.6);

header_row = [{'A'}, {'Ap'}];
write_data_cols = [paired_ctrst_scores, unpaired_ctrst_scores];

%writing data behind plot to file
xls_path = [paper_save_dir_sfig,  'contrast_scores_transitions_MBONA3.xls'];
[c] = write_xls_header(header_row, write_data_cols, xls_path);
write_data_cols = [];


%2. by odor identity
% PA_resps_paired = saved_resps_all2(1, :, 1)';
% PA_resps_unpaired = saved_resps_all2(1, :, 2)';
% PA_resps_ctrl = saved_resps_all2(1, :, 3)';
% BA_resps_paired = saved_resps_all2(2, :, 2)';
% BA_resps_unpaired = saved_resps_all2(2, :, 1)';
% BA_resps_ctrl = saved_resps_all2(2, :, 3)';
% 
% plot_mat2 = pad_n_concatenate(PA_resps_ctrl, PA_resps_paired, 2, nan);
% plot_mat2 = pad_n_concatenate(plot_mat2, PA_resps_unpaired, 2, nan);
% plot_mat2 = pad_n_concatenate(plot_mat2, BA_resps_ctrl, 2, nan);
% plot_mat2 = pad_n_concatenate(plot_mat2, BA_resps_paired, 2, nan);
% plot_mat2 = pad_n_concatenate(plot_mat2, BA_resps_unpaired, 2, nan);
% 
% paired_multiplier = 1;
% marker_colors = [PA_color; PA_color.*paired_multiplier; PA_color.*paired_multiplier;...
%     BA_color; BA_color.*paired_multiplier; BA_color.*paired_multiplier];
% line_colors = [];
% col_pairs = [];
% xlabels = [{'PA_c_t_r_l'}, {'PA_C_S_p_l'}, {'PA_C_S_m_n'}, {'BA_c_t_r_l'}, {'BA_C_S_p_l'}, {'BA_C_S_m_n'}];
% figure(8)
% [fig_h, r_vecs_saved] = scattered_dot_plot_ttest(plot_mat2, 8, 2.5, 4, 6.5, marker_colors, 1, col_pairs, line_colors, xlabels, 1, mean_color, 2, 0.05);
% ax_vals = axis;
% ax_vals(3) = -0.5;
% axis(ax_vals)
% ylabel('response size (dF/F)');
% fig_wrapup(fig_h, []);

% [hpairedPA, ppairedPA] = ttest(plot_mat2(:, 1), plot_mat2(:, 2))
% [hunpairedPA, punpairedPA] = ttest(plot_mat2(:, 1), plot_mat2(:, 3))
% [hpairedBA, ppairedBA] = ttest(plot_mat2(:, 4), plot_mat2(:, 5))
% [hunpairedBA, punpairedBA] = ttest(plot_mat2(:, 4), plot_mat2(:, 6))
