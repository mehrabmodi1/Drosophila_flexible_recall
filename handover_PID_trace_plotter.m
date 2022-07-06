clear all
close all

curr_dir = 'C:\Data\Data\Raw_data\20200203_PID\handover_PID_traces_set2\';

[del, odor_names1] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList.xls', 1);
[del, odor_names2] = xlsread('C:\Data\Code\general_code\IDnF_rig_code_20171031\Olfactometer\NewOlfactometer\calibration\odorList_olf2.xls', 1);
paper_save_dir_sfig = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\SFig2_PID_anemometer_traces\';


odor_names2{3} = 'Butyl acetate';
frame_time = 0.099;
[stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig] = load_params_trains_modular(curr_dir, [], frame_time);    %reading in trial stimulus parameters after matching time stamps to F traces

EL_color = [123,50,148]./256;
PA_color = [0,136,55]./256;
BA_color = [166,219,160]./256;


y_ax_lim = 0.04;

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

%plotting simple odor pulse PID traces
%PA
plot_simple_traces(3, 1, PA_color, 1, 2, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir, y_ax_lim, odor_names2, paper_save_dir_sfig)

%BA
plot_simple_traces(10, 3, BA_color, 3, 4, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir, y_ax_lim, odor_names2, paper_save_dir_sfig)

%EL
plot_simple_traces(11, 4, EL_color, 5, 6, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir, y_ax_lim, odor_names2, paper_save_dir_sfig)

%plotting handover odor pulse PID traces
%BA-BA
plot_hover_traces(10, 3, BA_color, BA_color, 13, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir, y_ax_lim, odor_names2, paper_save_dir_sfig)

%PA-BA
plot_hover_traces(3, 3, PA_color, BA_color, 7, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir, y_ax_lim, odor_names2, paper_save_dir_sfig)
%PA-EL
plot_hover_traces(3, 4, PA_color, EL_color, 8, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir, y_ax_lim, odor_names2, paper_save_dir_sfig)

%BA-PA
plot_hover_traces(10, 1, BA_color, PA_color, 9, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir, y_ax_lim, odor_names2, paper_save_dir_sfig)
%BA-EL
plot_hover_traces(10, 4, BA_color, EL_color, 10, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir, y_ax_lim, odor_names2, paper_save_dir_sfig)

%EL-PA
plot_hover_traces(11, 1, EL_color, PA_color, 11, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir, y_ax_lim, odor_names2, paper_save_dir_sfig)
%EL-BA
plot_hover_traces(11, 3, EL_color, BA_color, 12, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir, y_ax_lim, odor_names2, paper_save_dir_sfig)




function plot_simple_traces(olf1_odn, olf2_odn, od_color, fign1, fign2, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir, y_ax_lim, odor_names2, paper_save_dir_sfig)
frame_time = 0.099;

write_data_cols = [];

%list of all handover trials, with any combination of odors
handover_trs = find(stim_mat_simple(:, dur_col_ns(1)) == 10 & stim_mat_simple(:, dur_col_ns(2)) == 10);

%plotting olf1 traces
curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == olf1_odn & stim_mat_simple(:, dur_col_ns(1)) == 10);
[del1, del] = intersect(curr_trs, handover_trs);
curr_trs(del) = [];
curr_traces = get_PID_traces(curr_dir, curr_trs, frame_time, 0);
curr_traces(:, 1:2:end) = curr_traces(:, 1:2:end) - curr_traces(:, 2:2:end); 
curr_traces = curr_traces(1:(end - 5), 1:2:end);  %getting rid of LED traces
stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
stim_frs = [stim_frs{1}]; 
figure(fign1)
mean_tr = mean(curr_traces, 2, 'omitnan');
se_tr = std(curr_traces, [], 2, 'omitnan')./sqrt(size(curr_traces, 2));
%shadedErrorBar([], mean_tr, se_tr, {'Color', [0.6, 0.6, 0.6]}, 1);
plot(curr_traces, 'lineWidth', 1, 'Color', [0.6, 0.6, 0.6]);
hold on

ylabel('PID signal (V)');
ax_vals = axis;
ax_vals(4) = y_ax_lim;
axis(ax_vals);
set_xlabels_time(fign1, frame_time, 10);

n_traces1 = size(curr_traces, 2);
write_data_cols = pad_n_concatenate(write_data_cols, curr_traces, 2, nan);

%plotting olf2 traces
curr_trs = find(stim_mat_simple(:, od_col_ns(2)) == olf2_odn & stim_mat_simple(:, dur_col_ns(2)) == 10);
[del1, del] = intersect(curr_trs, handover_trs);
curr_trs(del) = [];
curr_trs(1) = [];
curr_traces = get_PID_traces(curr_dir, curr_trs, frame_time, 0);
curr_traces(:, 1:2:end) = curr_traces(:, 1:2:end) - curr_traces(:, 2:2:end); 
curr_traces = curr_traces(:, 1:2:end);  %getting rid of LED traces
stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
stim_frs = [stim_frs{2}]; 
mean_tr = mean(curr_traces, 2, 'omitnan');
se_tr = std(curr_traces, [], 2, 'omitnan')./sqrt(size(curr_traces, 2));
shadedErrorBar([], mean_tr, se_tr, {'Color', [0, 0, 0]}, 1);
plot(curr_traces, 'lineWidth', 1, 'Color', [0, 0, 0]);
fig_wrapup(fign1, [], [25, 30], 0.6);
add_stim_bar(fign1, stim_frs, od_color);

n_traces2 = size(curr_traces, 2);
write_data_cols = pad_n_concatenate(write_data_cols, curr_traces, 2, nan);

%writing data underlying plots to file
odname = odor_names2{olf2_odn};
col_header(1, 1:n_traces1) = {['OM1_', odname]};
col_header(1, (n_traces1 + 1):(n_traces1 + n_traces2)) = {['OM1_', odname]};
xls_path = [paper_save_dir_sfig,  'PID_traces_simple_', odname, '.xlsx'];
[c] = write_xls_header(col_header, write_data_cols, xls_path);
write_data_cols = [];
end


function [] = plot_hover_traces(olf1_odn, olf2_odn, od1_color, od2_color, fign, stim_mat, stim_mat_simple, od_col_ns, dur_col_ns, curr_dir, y_ax_lim , odor_names2, paper_save_dir_sfig)
frame_time = 0.099;

%list of all handover trials, with any combination of odors
handover_trs = find(stim_mat_simple(:, dur_col_ns(1)) == 10 & stim_mat_simple(:, dur_col_ns(2)) == 10);

curr_trs = find(stim_mat_simple(:, od_col_ns(1)) == olf1_odn & stim_mat_simple(:, dur_col_ns(1)) == 10 ...
    & stim_mat_simple(:, od_col_ns(2)) == olf2_odn & stim_mat_simple(:, dur_col_ns(2)) == 10);
curr_trs(1) = [];
curr_traces = get_PID_traces(curr_dir, curr_trs, frame_time, 0);
curr_traces(:, 1:2:end) = curr_traces(:, 1:2:end) - curr_traces(:, 2:2:end); 
curr_traces = curr_traces(:, 1:2:end);  %getting rid of LED traces
stim_frs = compute_stim_frs_modular(stim_mat, curr_trs(1), frame_time);
stim_frs = [stim_frs{1}; stim_frs{2}];
mean_tr = mean(curr_traces, 2, 'omitnan');
se_tr = std(curr_traces, [], 2, 'omitnan')./sqrt(size(curr_traces, 2));

figure(fign)
%shadedErrorBar([], mean_tr, se_tr, {'Color', [0.6, 0.6, 0.6]}, 1);
plot(curr_traces, 'lineWidth', 1, 'Color', [0.65, 0.65, 0.65]);
ylabel('PID signal (V)');
ax_vals = axis;
ax_vals(4) = y_ax_lim;
axis(ax_vals);
set_xlabels_time(fign, frame_time, 10);
fig_wrapup(fign, [], [25, 30], 0.6);
add_stim_bar(fign, stim_frs, [od1_color; od2_color]);

%only writing PABA or BAPA to disk
write_data_cols = curr_traces;
n_traces = size(curr_traces, 2);
if olf1_odn == 3 && olf2_odn == 3 
    %writing data underlying plots to file
    col_header(1, 1:n_traces) = {['OM1_PA_OM2_BA']};    
    xls_path = [paper_save_dir_sfig,  'PID_traces_transition_PABA.xlsx'];
    [c] = write_xls_header(col_header, write_data_cols, xls_path);
elseif olf1_odn == 10 && olf2_odn == 1
    col_header(1, 1:n_traces) = {['OM1_BA_OM2_PA']};
    xls_path = [paper_save_dir_sfig,  'PID_traces_transition_BAPA.xlsx'];
    [c] = write_xls_header(col_header, write_data_cols, xls_path);    
end

end


