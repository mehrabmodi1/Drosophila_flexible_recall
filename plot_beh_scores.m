clear all
close all

A3_path = 'C:\Data\Data\beh_data_extr_PIs\Alpha3_fine_coarse.xls';   %US-subst 10x training, 1 day memory 
A3_genpath = 'C:\Data\Data\beh_data_extr_PIs\Alpha3_generalization.xls';   %US-subst 10x training, 1 day memory 
MB296B_genpath = 'C:\Data\Data\beh_data_extr_PIs\MB296B_gen_data_3xtraining.xls';
MB296B_discrpath = 'C:\Data\Data\beh_data_extr_PIs\MB296B_discr_data_3xtraining.xls';
TH_rescue_discrpath = 'C:\Data\Data\beh_data_extr_PIs\THRescue_fine_discr.xls';


paper_save_dir = 'C:\Backup\Stuff\Janelia\paper_drafts\Mehrab_papers\PaBaEl2\fig_data\Fig1_behavior\';

%Plotting fine, coarse discr data for Alpha3, 3x training, 1 day memory
[alpha3_data] = xlsread(A3_path, 1).*-1;
marker_colors = repmat([0.65, 0.65, 0.65], 2, 1);
line_colors = [];
col_pairs = [];
%Note: this data is for these drivers: MB630B(B1B2), MB630B,  
xlabels = [{'A3PAEL'}, {'PABA'}] ;

figure(3)
fig_h = scattered_dot_plot_ttest(alpha3_data, 3, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
hold on
ax_vals = axis;
ax_vals(3) = -0.2;
ax_vals(4) = 0.85;
axis(ax_vals);
plot([0, ax_vals(2)], [0, 0], 'Color', [0. 0. 0]);
ylabel('1 day memory PI');
fig_wrapup(fig_h, [], [100, 25], 0.6);

%writing data behind plot to file
write_data_cols = alpha3_data;
header_row = [{'AvsB'}, {'ApvsB'}];
xls_path = [paper_save_dir,  'easy_hard_discr_Alpha3_MB630B.xls'];
[c] = write_xls_header(header_row, write_data_cols, xls_path);
write_data_cols = [];


%Plotting MB296B discr and gen scores separately
[discr_data_MB296B] = xlsread(MB296B_discrpath, 1);
figure(5)
xlabels = [{'G2AvsB'}, {'AvsAp'}];
fig_h = scattered_dot_plot_ttest(discr_data_MB296B, 5, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
hold on
ylabel('discrimination score (PI)');
ax_vals = axis;
ax_vals(3) = -0.2;
ax_vals(4) = 0.85;
axis(ax_vals);
plot([0, ax_vals(2)], [0, 0], 'Color', [0, 0, 0]);
fig_wrapup(fig_h, [], [100, 25], 0.6);

%writing data behind plot to file
write_data_cols = discr_data_MB296B;
header_row = [{'AvsB'}, {'AvsAp'}];
xls_path = [paper_save_dir,  'easy_hard_discr_G2Ap1_MB296B.xls'];
[c] = write_xls_header(header_row, write_data_cols, xls_path);
write_data_cols = [];



figure(6)
[gen_data_MB296B] = xlsread(MB296B_genpath, 1);
xlabels = [{'G2PAEL'}, {'BAEL'}];
fig_h = scattered_dot_plot_ttest(gen_data_MB296B, 6, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
hold on
plot([0, ax_vals(2)], [0, 0], 'Color', [0, 0, 0]);
ylabel('generalization score');
ax_vals = axis;
ax_vals(4) = 0.85;
ax_vals(1) = 0.8;
ax_vals(3) = -0.6;
axis(ax_vals);
fig_wrapup(fig_h, [], [100, 25], 0.6);

%writing data behind plot to file
write_data_cols = gen_data_MB296B;
header_row = [{'AvsB'}, {'ApvsB'}];
xls_path = [paper_save_dir,  'easydiscr_gen_G2Ap1_MB296B.xls'];
[c] = write_xls_header(header_row, write_data_cols, xls_path);
write_data_cols = [];


[gen_data_Alpha3] = xlsread(A3_genpath, 1);
figure(7)
xlabels = [{'A3PAEL'}, {'BAEL'}];
fig_h = scattered_dot_plot_ttest(gen_data_Alpha3, 7, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
hold on
plot([0, ax_vals(2)], [0, 0], 'Color', [0, 0, 0]);
ylabel('generalization score');
ax_vals = axis;
ax_vals(4) = 0.85;
ax_vals(1) = 0.8;
ax_vals(3) = -0.6;
axis(ax_vals);
fig_wrapup(fig_h, [], [100, 25], 0.6);

%writing data behind plot to file
write_data_cols = gen_data_Alpha3;
header_row = [{'AvsB'}, {'ApvsB'}];
xls_path = [paper_save_dir,  'easydiscr_gen_Alpha3_MB630B.xls'];
[c] = write_xls_header(header_row, write_data_cols, xls_path);
write_data_cols = [];


%Plotting TH-rescue, fine discr for MB296B data (only G2A'1)
[rescue_data] = xlsread(TH_rescue_discrpath, 1);
marker_colors = repmat([0.65, 0.65, 0.65], 2, 1);
line_colors = [];
col_pairs = [];
xlabels = [{'G2res'}, {'dblknk'}];
figure(8)
fig_h = scattered_dot_plot_ttest(rescue_data(:, [1:2]), 8, .6, 1, 4, marker_colors, 1, col_pairs, line_colors, xlabels, 2, mean_color, 2, 0.05, 0, 1, 'force_mean', wrapup_vars);
ylabel('performance index');
fig_wrapup(fig_h, [], [100, 25], 0.6);
ax_vals(3) = -0.4;
ax_vals(4) = 0.6;
ax_vals(1) = 0.8;
axis(ax_vals);
hold on
plot([0, ax_vals(2)], [0, 0], 'Color', [0, 0, 0]);
hold off

%writing data behind plot to file
write_data_cols = rescue_data(:, [1, 2]);
header_row = [{'AvsAp_rescue'}, {'AvsAp_norescue'}];
xls_path = [paper_save_dir,  'TH_rescue.xls'];
[c] = write_xls_header(header_row, write_data_cols, xls_path);
write_data_cols = [];

[p, h] = ranksum(rescue_data(:, 1), rescue_data(:, 2))


ax_vals = axis;
ax_vals(4) = 0.85;
axis(ax_vals);

%one-sample ttesting, checking if diff from 0
[p_MB296B_discr] = signrank_batch(discr_data_MB296B)
[p_MB296B_gen] = signrank_batch(gen_data_MB296B)
[p_screen] = signrank_batch(scr_data) 
[p_rescue] = signrank_batch(rescue_data(:, [1:2])) 
 p_alpha3 = signrank_batch(alpha3_data(:, [1:2])) 
 
 [p_alpha3_comparison, h] = ranksum(alpha3_data(:, 1), alpha3_data(:, 2));
 [p_alpha3] = bonf_holm([p_alpha3_comparison, p_alpha3])      %0.007 is the p value for the difference between easy and hard
 
 p_alpha3_MB296B_hard = ranksum(discr_data_MB296B(:, 2), alpha3_data(:, 2) );