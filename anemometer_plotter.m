clear all
close all

path = 'C:\Data\Data\Raw_data\20200803_anemometer\anemometer_acqn\';
paper_save_dir_sfig = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\SFig2_PID_anemometer_traces\';

frame_time = 0.099;
[traces, traces_orig] = get_PID_traces(path, [1:15], frame_time, 0);
[stim_mat, stim_mat_simple, column_heads, color_vec, good_tr_list, params_orig] = load_params_trains_modular(path, [], frame_time);

stim_latency = stim_mat(1).stimLatency;
acqn_time = abs(mean(diff(traces_orig(:, 3)), 'omitnan'));  %computing the PID acquisition rate from the time-stamp vector
fr_acq_ratio = frame_time./acqn_time;

%getting rid of LED acqn traces from traces and LED, timevec from traces_orig
traces = traces(:, 1:2:end);
traces_orig = traces_orig(:, 1:3:end);

%adding manually measured offset 
traces_orig = traces_orig + 0.155;  %155 mV  

%subtracting mean, pre-stimulus baseline from each trace
bk_vec = mean(traces_orig(1:round((stim_latency - 1)./acqn_time), :), 'omitnan');
bk_repmat = repmat(bk_vec, size(traces_orig, 1), 1);

%computing anemometer 'dF/F' or norm. traces
%traces_dbkbk = (traces_orig - bk_repmat)./bk_repmat;
traces_dbkbk = (traces_orig)./bk_repmat;
traces_f = movmean(traces_dbkbk, fr_acq_ratio, 1);

%Plotting
olf1_color = [0, 0, 0];
olf2_color = [0.65, 0.65, 0.65];


write_data_cols = [];
%olf1 trace
mean_trace = mean(traces_f(:, 2:5), 2, 'omitnan');
stim_frs = compute_stim_frs_modular(stim_mat, 2, acqn_time);
stim_frs = stim_frs{1} - 1000;
se_trace = std(traces_f(:, 2:5), [], 2, 'omitnan')./sqrt(4);
figure(1)
%shadedErrorBar([], mean_trace(1000:end), se_trace(1000:end), {'Color', [0, 0, 0]});
plot(traces_f(:, 2:5), 'Color', [0.65, 0.65, 0.65]);
hold on
ylabel('norm. air flow');
set_xlabels_time(1, acqn_time, 10);
ax_vals = axis;
plot([0, ax_vals(2)], [1, 1], '--', 'Color', [0.65, 0.65, 0.65])
ax_vals(3) = 0;
ax_vals(4) = 1.2;
axis(ax_vals);
fig_wrapup(1, [], [50, 60]);
od_color = olf1_color;
add_stim_bar(1, stim_frs, od_color);

write_data_cols = pad_n_concatenate(write_data_cols, mean(traces_f(:, 2:5), 2, 'omitnan'), 2, nan);


%Plotting
%olf2 trace
mean_trace = mean(traces_f(:, 6:10), 2, 'omitnan');
stim_frs = compute_stim_frs_modular(stim_mat, 6, acqn_time);
stim_frs = stim_frs{2} - 1000;
se_trace = std(traces_f(:, 6:10), [], 2, 'omitnan')./sqrt(4);
figure(2)
%shadedErrorBar([], mean_trace(1000:end), se_trace(1000:end), {'Color', [0, 0, 0]});
plot(traces_f(:, 6:10), 'Color', [0.65, 0.65, 0.65]);
ylabel('norm. air flow');
set_xlabels_time(2, acqn_time, 10);
axis(ax_vals);
hold on
plot([0, ax_vals(2)], [1, 1], '--', 'Color', [0.65, 0.65, 0.65])
fig_wrapup(2, [], [50, 60]);
od_color = olf2_color;
add_stim_bar(2, stim_frs, od_color);

write_data_cols = pad_n_concatenate(write_data_cols, mean(traces_f(:, 6:10), 2, 'omitnan'), 2, nan);


%Plotting
%olf1-olf2 handover trace
mean_trace = mean(traces_f(:, 11:15), 2, 'omitnan');
stim_frs = compute_stim_frs_modular(stim_mat, 11, acqn_time);
stim_frs = [stim_frs{1}; stim_frs{2}] - 1000;
se_trace = std(traces_f(:, 11:15), [], 2, 'omitnan')./sqrt(4);
figure(3)
%shadedErrorBar([], mean_trace(1000:end), se_trace(1000:end), {'Color', [0, 0, 0]});
plot(traces_f(:, 11:15), 'Color', [0.65, 0.65, 0.65]);
ylabel('norm. air flow');
set_xlabels_time(3, acqn_time, 10);
ax_vals(2) = 100000;
hold on
plot([0, ax_vals(2)], [1, 1], '--', 'Color', [0.65, 0.65, 0.65])
axis(ax_vals);
fig_wrapup(3, [], [50, 60]);
od_color = [olf1_color; olf2_color];
add_stim_bar(3, stim_frs, od_color);

write_data_cols = pad_n_concatenate(write_data_cols, mean(traces_f(:, 11:15), 2, 'omitnan'), 2, nan);



%writing data underlying plots to file
col_heads = [{'odormachine1'}, {'odormachine2'}, {'OM1toOM2'}];
xls_path = [paper_save_dir_sfig,  'anemometer_traces.xlsx'];
[c] = write_xls_header(col_heads, write_data_cols, xls_path);
write_data_cols = [];

