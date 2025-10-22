function results = analyze_RSTS(varargin)
% =========================================================================
%                COMPREHENSIVE RAMAN & THOMSON ANALYSIS (V9.9 GUI - RSTSGUI2)
% =========================================================================
% V9.9 Changes (RSTSGUI2):
% - ENHANCED: Weighted fitting in Int.Calib.Full prioritizes central region
%   around 532 nm using Gaussian weights (sigma=5 nm) to prevent fits from
%   being pulled by distant wings
% - ENHANCED: Monte Carlo-based error estimation applied to ALL three fitting
%   methods (Int.Calib.Full, Shape Calib., Int.Calib.Gauss) provides realistic
%   parameter uncertainties by finding all Te/ne combinations that produce
%   visually acceptable fits (within 1.5x best SSR threshold)
% - IMPROVED: Errors now typically 10-100x larger, matching empirical ±15%
%   criterion for visually indistinguishable synthetic spectra
% - UNIFIED: All fitting callbacks now use consistent Monte Carlo error approach
%
% V9.8 Changes:
% - ADDED: Automatic update of Te, ne, alpha fields after Gaussian fit
% - ADDED: Gamma field is now editable
% - FIXED: "Recalculate TS" button functionality simplified
% - ADDED: Status message display in GUI
% - ENHANCED: "Recalibrate from Raman" uses full multi-start optimization
% - ADDED: Number of multi-starts control for Raman and Shape optimization
% - ADDED: Save buttons for each axis to save plots as JPG
%
% V9.7 Changes:
% - FIXED: Added the missing helper function 'coherent_intensity_model' to
%   enable the "Optimize All" iterative fit and resolve the "Undefined
%   function" error.
% - FIXED: Added the missing 'calculate_alpha_from_ne' helper function to
%   enable the "Recalculate TS" button functionality.
% =========================================================================

%% --- 1. CHOOSE DEFAULT FILE SET ---
% 1: NYU Data
% 2: CC Plasma-Liquid
% 3: E-Beam (Texas A&M)
% 4: RF Jet (WUSTL)
default_set = 2;

%% --- 2. DEFINE FILE PATHS AND PARAMETERS BASED ON FLAG ---
switch default_set
    case 1 % NYU Data
        disp('Loading Default Set 1: NYU Data');
        % p.raman_file = 'G:\My Drive\data\NYU\09.02.2025 TS\250 us Pos Pol\gr2400 sw150 um PolHorz\Raman\Raman 300 gpe AvImage.tif';
        % p.thomson_file = 'G:\My Drive\data\NYU\09.02.2025 TS\250 us Pos Pol\gr2400 sw150 um PolHorz\5 mA\afternoon\vary h2l\hpl 15AvImage.tif';
        % p.raman_bg_file = ''; p.thomson_bg_file = '';
        % p.yc = 90; p.wy = 40; p.exclusion_cols = 455:530;
        % p.shots_raman = 300; p.shots_thomson = 300;
        % p.energy_J_raman = 52e-3; p.energy_J_thomson = 52e-3;
        % p.pressure_Pa = 101325; p.temp_K_initial = 293.15;
        % p.a_initial = 0.010050; p.b_initial = 526.7384;
        % p.linewidth_initial = 0.1004; p.mixing_ratios = [0.79, 0.21];

        p.raman_file = 'G:\My Drive\data\NYU\09.12.2025 TS 2.5 mS\Raman\Raman_AvImage.tif';
        p.thomson_file = 'G:\My Drive\data\NYU\09.12.2025 TS 2.5 mS\5 mA\TS\RAW\-15 h2l 65 p2l 44 p2w 400 gpe gain5K_AvImage.tif';
        p.raman_bg_file = ''; p.thomson_bg_file = '';
        p.yc = 80; p.wy = 40; p.exclusion_cols = 455:550;
        p.shots_raman = 400; p.shots_thomson = 400;
        p.energy_J_raman = 52e-3; p.energy_J_thomson = 52e-3;
        p.pressure_Pa = 101325; p.temp_K_initial = 293.15;
        p.a_initial = 0.010050; p.b_initial = 526.7384;
        p.linewidth_initial = 0.1004; p.mixing_ratios = [0.79, 0.21];


    case 2 % CC Plasma-Liquid
        disp('Loading Default Set 2: CC Plasma-Liquid Data');
        folder = 'G:\My Drive\data\CC plasma-liquid\10.31.2024 LTS 500 us cond\t ev\after lunch\366140 ns\';
        p.thomson_file = [folder, '366220 ns 7.85 mm l2p 4.04 h2l  ChB 366 us 532 nm gd 100 ns gw 30 ns Gain 10K gpe 3000 14_51_48.tif'];
        RSfolder='G:\My Drive\data\CC plasma-liquid\10.31.2024 LTS 500 us cond\t ev\after lunch\Raman\';
        p.raman_file   = [RSfolder, 'Raman 8.45 mm l2p 4.04 h2l  ChB 366 us 532 nm gd 100 ns gw 30 ns Gain 10K gpe 3000 16_04_49.tif'];
        bgfolder='G:\My Drive\data\CC plasma-liquid\10.31.2024 LTS 500 us cond\t ev\after lunch\bg\';
        bgfile='bg 8.45 mm l2p 4.04 h2l  ChB 366 us 532 nm gd 100 ns gw 30 ns Gain 10K gpe 3000 16_02_04.tif';
        p.raman_bg_file = [bgfolder, bgfile]; p.thomson_bg_file = [bgfolder, bgfile];
        p.yc = 46; p.wy = 40; p.exclusion_cols = 475:490;
        p.shots_raman = 3000; p.shots_thomson = 3000;
        p.energy_J_raman = 12e-3; p.energy_J_thomson = 12e-3;
        p.pressure_Pa = 101325; p.temp_K_initial = 293.15;
        p.a_initial = 0.01311544401544386; p.b_initial = 525.702629150;
        p.linewidth_initial = 0.067; p.mixing_ratios = [0.79, 0.21];

    case 3 % E-Beam Data
        disp('Loading Default Set 3: E-Beam Data (Texas A&M)');
        TSfolder = 'G:\My Drive\data\TexasA&M Junhwi Bak Feb2025\02.13.2025\Argon\';
        p.thomson_file = [TSfolder, 'lts_110G_acc1000_181mJ_1p46mTorr 15_02_26 1.tiff'];
        RSfolder='G:\My Drive\data\TexasA&M Junhwi Bak Feb2025\02.13.2025\RRS_Calib_after\';
        p.raman_file   = [RSfolder, 'RRS_202mTorr_171mJ_acc1000 17_44_09 1.tiff'];
        bgfolderRS='G:\My Drive\data\TexasA&M Junhwi Bak Feb2025\02.13.2025\RRS_Calib_after\';
        bgfileRS='bg_202mTorr_171mJ_acc1000 17_45_10 1.tiff';
        bgfolderTS='G:\My Drive\data\TexasA&M Junhwi Bak Feb2025\02.13.2025\Argon\';
        bgfileTS='pla_110G_acc1000_181mJ_1p46mTorr 15_01_23 1.tiff';
        p.raman_bg_file = [bgfolderRS, bgfileRS]; p.thomson_bg_file = [bgfolderTS, bgfileTS];
        p.yc = 700; p.wy = 230; p.exclusion_cols = 524:596;
        p.shots_raman = 1000; p.shots_thomson = 1000;
        p.energy_J_raman = 171e-3; p.energy_J_thomson = 181e-3;
        p.pressure_Pa = 26931.1; p.temp_K_initial = 293.15;
        p.a_initial = 0.029937; p.b_initial = 515.2500;
        p.linewidth_initial = 0.7409; p.mixing_ratios = [1, 0];

    case 4 % RF Jet WUSTL Data
        disp('Loading Default Set 4: RF Jet Data (WUSTL)');
        TSfolder = 'G:\My Drive\data\RF Jet WUSTL\11.02.20 thomson hal 8 mm\thomson\nov 2\';
        p.thomson_file = [TSfolder, 'thomson le 64 mj  10000 shots 3088 2020 November 02 17_29_44.tif'];
        RSfolder='G:\My Drive\data\RF Jet WUSTL\09.30.2020 Bg and RAMAN for 10000 shots\Raman\';
        p.raman_file   = [RSfolder, 'raman lp 80 mj 10000 shots 2020 September 30 16_56_32.tif'];
        bgfolderRS='G:\My Drive\data\RF Jet WUSTL\09.30.2020 Bg and RAMAN for 10000 shots\bg  plasma on laser off\';
        bgfileRS='bg l-off p-on 10000 2020 September 30 17_17_56.tif';
        bgfolderTS='G:\My Drive\data\RF Jet WUSTL\09.30.2020 Bg and RAMAN for 10000 shots\bg  plasma on laser off\';
        bgfileTS='bg l-off p-on 10000 2020 September 30 17_17_56.tif';
        p.raman_bg_file = [bgfolderRS, bgfileRS]; p.thomson_bg_file = [bgfolderTS, bgfileTS];
        p.yc = 532; p.wy = 60; p.exclusion_cols = 506:529;
        p.shots_raman = 10000; p.shots_thomson = 10000;
        p.energy_J_raman = 80e-3; p.energy_J_thomson = 64e-3;
        p.pressure_Pa = 101325; p.temp_K_initial = 295.2;
        p.a_initial = 0.0124919; p.b_initial = 525.5345;
        p.linewidth_initial = 0.0669; p.mixing_ratios = [0.79, 0.21];

    otherwise
        error('Invalid default_set flag. Please choose 1, 2, 3, or 4.');
end

%% --- 3. PARSE INPUT ARGUMENTS, LOAD DATA, AND VALIDATE PARAMS ---
parser = inputParser;
addParameter(parser, 'raman_file', p.raman_file, @ischar);
addParameter(parser, 'thomson_file', p.thomson_file, @ischar);
addParameter(parser, 'raman_bg_file', p.raman_bg_file, @ischar);
addParameter(parser, 'thomson_bg_file', p.thomson_bg_file, @ischar);
addParameter(parser, 'yc', p.yc, @isnumeric);
addParameter(parser, 'wy', p.wy, @isnumeric);
addParameter(parser, 'exclusion_cols', p.exclusion_cols, @isnumeric);
addParameter(parser, 'shots_raman', p.shots_raman, @isnumeric);
addParameter(parser, 'shots_thomson', p.shots_thomson, @isnumeric);
addParameter(parser, 'energy_J_raman', p.energy_J_raman, @isnumeric);
addParameter(parser, 'energy_J_thomson', p.energy_J_thomson, @isnumeric);
addParameter(parser, 'laser_duration_s', 7e-9, @isnumeric);
addParameter(parser, 'pressure_Pa', p.pressure_Pa, @isnumeric);
addParameter(parser, 'temp_K_initial', p.temp_K_initial, @isnumeric);
addParameter(parser, 'a_initial', p.a_initial, @isnumeric);
addParameter(parser, 'b_initial', p.b_initial, @isnumeric);
addParameter(parser, 'linewidth_initial', p.linewidth_initial, @isnumeric);
addParameter(parser, 'mixing_ratios', p.mixing_ratios, @isnumeric);

parse(parser, varargin{:});
params = parser.Results;

params.raman_img_raw = double(imread(params.raman_file));
if ~isempty(params.raman_bg_file), params.raman_img_raw = params.raman_img_raw - double(imread(params.raman_bg_file)); end
params.thomson_img_raw = double(imread(params.thomson_file));
if ~isempty(params.thomson_bg_file), params.thomson_img_raw = params.thomson_img_raw - double(imread(params.thomson_bg_file)); end

[nRows, nCols] = size(params.raman_img_raw);
params_adjusted = false;
yc = params.yc; wy = params.wy;
if (yc + floor(wy/2)) > nRows || (yc - floor(wy/2)) < 1
    params.yc = round(nRows / 2);
    params.wy = min(wy, floor(nRows/2));
    if mod(params.wy, 2) ~= 0, params.wy = params.wy - 1; end
    params_adjusted = true;
end
excl_start = min(params.exclusion_cols); excl_end = max(params.exclusion_cols);
if excl_start < 1 || excl_end > nCols
    new_center_col = round(nCols / 2);
    new_excl_width = min(excl_end - excl_start + 1, round(nCols / 10));
    new_start = new_center_col - floor(new_excl_width / 2);
    new_end = new_center_col + floor(new_excl_width / 2);
    params.exclusion_cols = new_start:new_end;
    params_adjusted = true;
end
if params_adjusted
    fprintf('INFO: Default ROI/Exclusion parameters were auto-adjusted to fit initial image dimensions.\n');
end

params.laser_wavelength = 532.0;

%% --- 4. INITIAL ANALYSIS ---
fprintf('=== Performing Initial Analysis ===\n');
params.signal_rows = (params.yc - floor(params.wy/2)):(params.yc + floor(params.wy/2));
[results, plot_data] = run_full_analysis(params);

%% --- 5. LAUNCH MAIN GUI ---
fprintf('\nLaunching main analysis GUI...\n');
generate_summary_gui(plot_data, params);

fprintf('\nAnalysis complete. Close figures to exit.\n');
end

% =========================================================================
%                        MAIN ANALYSIS PIPELINE
% =========================================================================
function [results, plot_data] = run_full_analysis(params)
me_c2_eV = 0.511e6;
sigma_T_differential = 7.94e-30;

fprintf('\n--- Analyzing Raman Data ---\n');
raman_img_processed = block_central_pixels(params.raman_img_raw, params.signal_rows, params.exclusion_cols);
pixels_exp = 1:size(params.raman_img_raw, 2);
experimental_raman_spectrum = mean(raman_img_processed(params.signal_rows, :), 1);
experimental_raman_spectrum = subtract_wing_baseline(experimental_raman_spectrum);

species_list = {'N2', 'O2'};
mixing_ratios = params.mixing_ratios;

wavelength_exp = params.a_initial * pixels_exp + params.b_initial;
fitted_linewidth = params.linewidth_initial;
theoretical_spectrum = calculate_theoretical_spectrum(wavelength_exp, species_list, mixing_ratios, params.pressure_Pa, params.temp_K_initial, fitted_linewidth, params.laser_wavelength);
experimental_integral = trapz(wavelength_exp, experimental_raman_spectrum);
theoretical_integral = trapz(wavelength_exp, theoretical_spectrum);
total_energy_raman = params.energy_J_raman * params.shots_raman;
if theoretical_integral > 0
    calibration_factor = experimental_integral / (total_energy_raman * theoretical_integral);
else
    calibration_factor = 1e-15;
end

raman_results.wavelength_axis = wavelength_exp;
raman_results.experimental_spectrum = experimental_raman_spectrum;
raman_results.theoretical_spectrum_model = theoretical_spectrum;
raman_results.calibration_factor = calibration_factor;

fprintf('--- Analyzing Thomson Data ---\n');
thomson_img_processed = block_central_pixels(params.thomson_img_raw, params.signal_rows, params.exclusion_cols);
raw_thomson_spectrum = mean(thomson_img_processed(params.signal_rows, :), 1);
wavelength_nm = params.a_initial * pixels_exp + params.b_initial;
thomson_spectrum_processed = subtract_wing_baseline(raw_thomson_spectrum);

% Multistart fit
fit_mask = ~ismember(1:length(pixels_exp), params.exclusion_cols);
gauss_model = @(p, x) p(1) * exp(-4 * log(2) * ((x - p(2)).^2 / p(3)^2)) + p(4);

% Multi-start Gaussian fitting
if isfield(params, 'num_starts') && ~isempty(params.num_starts)
    num_starts = params.num_starts;
else
    num_starts = 6; % Fallback default for initial analysis
end
best_p_fit_gauss = [];
best_resnorm_gauss = inf;

for start_idx = 1:num_starts
    if start_idx == 1
        p0_gauss = [max(thomson_spectrum_processed), params.laser_wavelength, 2.0, 0];
    else
        p0_gauss = [max(thomson_spectrum_processed)*(0.8+0.4*rand()), params.laser_wavelength+0.2*randn(), 1.0+2.0*rand(), 0];
    end

    [p_fit_current, resnorm_current] = lsqcurvefit(gauss_model, p0_gauss, wavelength_nm(fit_mask), thomson_spectrum_processed(fit_mask), [], [], optimset('Display','final'));

    if resnorm_current < best_resnorm_gauss
        best_p_fit_gauss = p_fit_current;
        best_resnorm_gauss = resnorm_current;
    end
end

p_fit_gauss = best_p_fit_gauss;
resnorm_gauss = best_resnorm_gauss;

fwhm_uncorrected = p_fit_gauss(3); amp_uncorrected = p_fit_gauss(1);
Te_uncorrected = me_c2_eV / (32 * log(2)) * (fwhm_uncorrected / params.laser_wavelength)^2;
total_energy_thomson = params.energy_J_thomson * params.shots_thomson;
gauss_fit_spectrum = gauss_model(p_fit_gauss, wavelength_nm(fit_mask));
total_counts_uncorrected = sum(gauss_fit_spectrum);
ne_uncorrected = total_counts_uncorrected / (calibration_factor * total_energy_thomson * sigma_T_differential);

r_squared_gauss = 1 - (resnorm_gauss / sum((thomson_spectrum_processed(fit_mask) - mean(thomson_spectrum_processed(fit_mask))).^2));
fprintf('Gaussian Fit: SSR = %.2e, R² = %.3f\n', resnorm_gauss, r_squared_gauss);

results.uncorrected.Te_eV = Te_uncorrected;
results.uncorrected.ne_m3 = ne_uncorrected;

deconv_results = amplitude_fit_deconvolution(wavelength_nm, thomson_spectrum_processed, ...
    calibration_factor, fitted_linewidth, params, p_fit_gauss(3), p_fit_gauss(1));

results.corrected.Te_eV = deconv_results.Te_eV_corrected;
results.corrected.ne_m3 = deconv_results.ne_corrected;
results.calibration_factor = calibration_factor;
results.deconv_results = deconv_results;

plot_data = struct('thomson_img_processed', thomson_img_processed, 'raman_img_processed', raman_img_processed, ...
    'raman_results', raman_results, 'wavelength_nm', wavelength_nm, 'mfp', 8, ...
    'thomson_spectrum_processed', thomson_spectrum_processed, 'gauss_fit_params', p_fit_gauss, ...
    'uncorrected_results', results.uncorrected, 'r_squared_gauss', r_squared_gauss, ...
    'ssr_gauss', resnorm_gauss, ...
    'deconv_results', deconv_results);
end

% =========================================================================
%                        LOCAL PLOTTING & GUI FUNCTIONS
% =========================================================================

function generate_summary_gui(data, params)

% Adaptive sizing code
screen_size = get(0, 'ScreenSize');
screen_width = screen_size(3);
screen_height = screen_size(4);

% Calculate usable screen height (accounting for taskbar ~40-80px)
usable_height = screen_height - 80;

% Use 95% of screen size for better visibility
fig_width = min(screen_width * 0.95, 1600);
fig_height = min(usable_height * 0.95, 1100);

% Center the figure on screen
fig_left = max((screen_width - fig_width) / 2, 10);
fig_bottom = max((usable_height - fig_height) / 2, 10);

fig = uifigure('Name', 'RSTS Analysis V9.8', ...
    'Position', [fig_left, fig_bottom, fig_width, fig_height], ...
    'NumberTitle', 'off', ...
    'AutoResizeChildren', 'off', ...
    'Scrollable', 'on');

% Make the figure resizable
fig.Resize = 'on';
fig.WindowState = 'normal';  % Ensure it's not maximized initially

handles = struct();

handles.params = params;
handles.plot_data = data;
handles.current_ts_spectrum = data.thomson_spectrum_processed;
handles.current_results.uncorrected = data.uncorrected_results;
handles.current_deconv = data.deconv_results;

handles.raman_file_original = params.raman_file;
handles.thomson_file_original = params.thomson_file;
handles.raman_bg_file = params.raman_bg_file;
handles.thomson_bg_file = params.thomson_bg_file;
handles.plot_style = 'line';
handles.fit_pixel_shift = 0;

handles.temp_K_current = params.temp_K_initial;
handles.n2_fraction = params.mixing_ratios(1);
handles.o2_fraction = params.mixing_ratios(2);

% Main grid: plots on top, controls on bottom
main_grid = uigridlayout(fig, [2, 1]);
main_grid.RowHeight = {'3x', '1x'};  % More space for plots
main_grid.Scrollable = 'on';
main_grid.ColumnWidth = {'1x'};

% Make figure scrollable with explicit size
fig.Scrollable = 'on';
fig.HandleVisibility = 'on';
% TOP SECTION: Plots and Thomson Files
plot_panel = uipanel(main_grid, 'BorderType', 'none');
plot_panel.Layout.Row = 1;

% Create grid for plot panel: plots on left, file browser on right
plot_outer_grid = uigridlayout(plot_panel, [2, 2]);
plot_outer_grid.RowHeight = {'1x', 'fit'};  % Plots, then save buttons
plot_outer_grid.ColumnWidth = {'4x', '1x'};  % Plots wide, files narrow

% Left side: Create 2x2 grid for the 4 plots
plot_grid = uigridlayout(plot_outer_grid, [2, 2]);
plot_grid.Layout.Row = 1;
plot_grid.Layout.Column = 1;

% Right side: Thomson Files panel
file_browser_panel = uipanel(plot_outer_grid, 'Title', 'Thomson Files');
file_browser_panel.Layout.Row = 1;
file_browser_panel.Layout.Column = 2;

% Create file browser layout inside the panel
file_gl = uigridlayout(file_browser_panel);
file_gl.RowHeight = {'fit', '1x', 'fit'};
file_gl.ColumnWidth = {'1x'};

lbl = uilabel(file_gl, 'Text', 'Files in current folder:', 'FontWeight', 'bold');
lbl.Layout.Row = 1;
lbl.Layout.Column = 1;

handles.file_listbox = uilistbox(file_gl, 'ValueChangedFcn', @file_selected_callback, 'Multiselect', 'off');
handles.file_listbox.Layout.Row = 2;
handles.file_listbox.Layout.Column = 1;

btn = uibutton(file_gl, 'Text', 'Refresh File List', 'ButtonPushedFcn', @refresh_file_list_callback);
btn.Layout.Row = 3;
btn.Layout.Column = 1;

% Save buttons row (spans both columns at bottom)
save_button_grid = uigridlayout(plot_outer_grid, [1, 4]);
save_button_grid.Layout.Row = 2;
save_button_grid.Layout.Column = [1, 2];  % Span both columns
save_button_grid.ColumnWidth = repmat({'1x'}, 1, 4);

% Create the 4 axes in the plot_grid
handles.ax1 = uiaxes(plot_grid); 
handles.ax1.Layout.Row = 1; handles.ax1.Layout.Column = 1;

handles.ax2 = uiaxes(plot_grid);
handles.ax2.Layout.Row = 1; handles.ax2.Layout.Column = 2;

handles.ax3 = uiaxes(plot_grid); 
handles.ax3.Layout.Row = 2; handles.ax3.Layout.Column = 1;

handles.ax4 = uiaxes(plot_grid);
handles.ax4.Layout.Row = 2; handles.ax4.Layout.Column = 2;

% Create save buttons
handles.save_btn1 = uibutton(save_button_grid, 'Text', 'Save Raman Image', 'ButtonPushedFcn', @(src,evt)save_axis_plot(handles.ax1, 'raman_image'));
handles.save_btn2 = uibutton(save_button_grid, 'Text', 'Save Thomson Image', 'ButtonPushedFcn', @(src,evt)save_axis_plot(handles.ax2, 'thomson_image'));
handles.save_btn3 = uibutton(save_button_grid, 'Text', 'Save Raman Spectrum', 'ButtonPushedFcn', @(src,evt)save_axis_plot(handles.ax3, 'raman_spectrum'));
handles.save_btn4 = uibutton(save_button_grid, 'Text', 'Save TS Spectrum', 'ButtonPushedFcn', @(src,evt)save_axis_plot(handles.ax4, 'ts_spectrum'));

% Plot the initial data
handles.h_raman_img = imagesc(handles.ax1, data.raman_img_processed); 
colormap(handles.ax1, 'jet'); hold(handles.ax1, 'on');
handles.h_raman_roi_rect = rectangle(handles.ax1, 'Position', [1, params.signal_rows(1), size(data.raman_img_processed, 2)-1, length(params.signal_rows)], 'EdgeColor', 'y', 'LineWidth', 2, 'LineStyle', '--');
handles.h_raman_excl_rect = rectangle(handles.ax1, 'Position', [min(params.exclusion_cols), 1, length(params.exclusion_cols), size(data.raman_img_processed, 1)], 'EdgeColor', 'r', 'LineWidth', 2, 'FaceColor', 'r', 'FaceAlpha', 0.3);
title(handles.ax1, 'Raman Image (Raw)'); xlabel(handles.ax1, 'Pixel Column'); colorbar(handles.ax1); axis(handles.ax1, 'image');

handles.h_thomson_img = imagesc(handles.ax2, medfilt2(data.thomson_img_processed, [data.mfp data.mfp]));
colormap(handles.ax2, 'jet'); hold(handles.ax2, 'on');
handles.h_thomson_roi_rect = rectangle(handles.ax2, 'Position', [1, params.signal_rows(1), size(data.thomson_img_processed, 2)-1, length(params.signal_rows)], 'EdgeColor', 'y', 'LineWidth', 2, 'LineStyle', '--');
handles.h_excl_rect = rectangle(handles.ax2, 'Position', [min(params.exclusion_cols), 1, length(params.exclusion_cols), size(data.thomson_img_processed, 1)], 'EdgeColor', 'r', 'LineWidth', 2, 'FaceColor', 'r', 'FaceAlpha', 0.3);
title(handles.ax2, 'Thomson Image (Raw)'); xlabel(handles.ax2, 'Pixel Column'); ylabel(handles.ax2, 'Pixel Row'); colorbar(handles.ax2); axis(handles.ax2, 'image');

handles.h_raman_spec = plot(handles.ax3, data.raman_results.wavelength_axis, data.raman_results.experimental_spectrum, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Experimental'); hold(handles.ax3, 'on');
model_scaled = data.raman_results.theoretical_spectrum_model * (max(data.raman_results.experimental_spectrum)/max(data.raman_results.theoretical_spectrum_model));
handles.h_raman_model = plot(handles.ax3, data.raman_results.wavelength_axis, model_scaled, 'r--', 'LineWidth', 2, 'DisplayName', 'Model (scaled)');
handles.h_calib_title = title(handles.ax3, sprintf('Raman Spectrum (Calib. Factor: %.2e)', data.raman_results.calibration_factor));
xlabel(handles.ax3, 'Wavelength (nm)'); ylabel(handles.ax3, 'Intensity (counts)'); legend(handles.ax3); grid(handles.ax3, 'on'); xlim(handles.ax3, 'auto');

gauss_model = @(p, x) p(1) * exp(-4 * log(2) * ((x - p(2)).^2 / p(3)^2)) + p(4);
fit_mask = ~ismember(1:length(data.wavelength_nm), params.exclusion_cols);
handles.h_data_points = plot(handles.ax4, data.wavelength_nm(fit_mask), data.thomson_spectrum_processed(fit_mask), 'b-', 'DisplayName', 'Data'); hold(handles.ax4, 'on');
handles.h_gauss_fit = plot(handles.ax4, data.wavelength_nm, gauss_model(data.gauss_fit_params, data.wavelength_nm), 'r-', 'LineWidth', 2, 'DisplayName', 'Gaussian Fit');
set(handles.h_gauss_fit, 'UserData', get(handles.h_gauss_fit, 'YData'));
handles.h_title = title(handles.ax4, sprintf('Direct TS Gaussian Fit (SSR=%.2e, R^2=%.3f)\nTe = %.2f \\pm %.2f eV, ne = %.2e \\pm %.2e m^{-3}', data.ssr_gauss, data.r_squared_gauss, data.deconv_results.Te_eV_corrected, data.deconv_results.Te_eV_error, data.deconv_results.ne_corrected, data.deconv_results.ne_error));
xlabel(handles.ax4, 'Wavelength (nm)'); ylabel(handles.ax4, 'Intensity (counts)'); legend(handles.ax4, 'Location', 'best'); grid(handles.ax4, 'on');

% BOTTOM SECTION: Control Panel
control_panel_outer = uipanel(main_grid, 'BorderType', 'none', 'Scrollable', 'on');
control_panel_outer.Layout.Row = 2;

% Control panel grid: main controls on left, utilities on right
outer_gl = uigridlayout(control_panel_outer, [1, 2]);
outer_gl.ColumnWidth = {'3x', '2x'};

main_controls_panel = uipanel(outer_gl, 'BorderType', 'none');
main_controls_panel.Layout.Column = 1;

utilities_panel = uipanel(outer_gl, 'Title', 'Utilities');
utilities_panel.Layout.Column = 2;

% Main controls grid layout - 7 rows, 9 columns
gl = uigridlayout(main_controls_panel);
gl.RowHeight = repmat({'fit'}, 1, 7);
gl.ColumnWidth = repmat({'1x'}, 1, 9);
gl.Padding = [2 2 2 2];
gl.RowSpacing = 2;
gl.ColumnSpacing = 3;

% ROW 1: ROI yc, wy, Excl Start/End, Refresh button
lbl = uilabel(gl, 'Text','ROI yc:'); lbl.Layout.Row=1; lbl.Layout.Column=1;
handles.edit_yc = uieditfield(gl, 'numeric', 'Value', params.yc); handles.edit_yc.Layout.Row=1; handles.edit_yc.Layout.Column=2;
lbl = uilabel(gl, 'Text','wy:'); lbl.Layout.Row=1; lbl.Layout.Column=3;
handles.edit_wy = uieditfield(gl, 'numeric', 'Value', params.wy); handles.edit_wy.Layout.Row=1; handles.edit_wy.Layout.Column=4;
lbl = uilabel(gl, 'Text','Excl:'); lbl.Layout.Row=1; lbl.Layout.Column=5;
handles.edit_start = uieditfield(gl, 'numeric', 'Value', min(params.exclusion_cols)); handles.edit_start.Layout.Row=1; handles.edit_start.Layout.Column=6;
lbl = uilabel(gl, 'Text','to'); lbl.Layout.Row=1; lbl.Layout.Column=7;
handles.edit_end = uieditfield(gl, 'numeric', 'Value', max(params.exclusion_cols)); handles.edit_end.Layout.Row=1; handles.edit_end.Layout.Column=8;
btn = uibutton(gl, 'Text', 'Refresh', 'ButtonPushedFcn', @refresh_images_roi_callback, 'FontWeight', 'bold', 'BackgroundColor', [1, 0.95, 0.8]);
btn.Layout.Row=1; btn.Layout.Column=9;

% ROW 2: Shots, Coefficients, Toggle
lbl = uilabel(gl, 'Text','RS Shots:'); lbl.Layout.Row=2; lbl.Layout.Column=1;
handles.edit_shots_r = uieditfield(gl, 'numeric', 'Value', params.shots_raman); handles.edit_shots_r.Layout.Row=2; handles.edit_shots_r.Layout.Column=2;
lbl = uilabel(gl, 'Text','TS Shots:'); lbl.Layout.Row=2; lbl.Layout.Column=3;
handles.edit_shots_t = uieditfield(gl, 'numeric', 'Value', params.shots_thomson); handles.edit_shots_t.Layout.Row=2; handles.edit_shots_t.Layout.Column=4;
lbl = uilabel(gl, 'Text','a:'); lbl.Layout.Row=2; lbl.Layout.Column=5;
handles.edit_a = uieditfield(gl, 'numeric', 'Value', params.a_initial, 'ValueDisplayFormat', '%.6f'); handles.edit_a.Layout.Row=2; handles.edit_a.Layout.Column=6;
lbl = uilabel(gl, 'Text','b:'); lbl.Layout.Row=2; lbl.Layout.Column=7;
handles.edit_b = uieditfield(gl, 'numeric', 'Value', params.b_initial, 'ValueDisplayFormat', '%.4f'); handles.edit_b.Layout.Row=2; handles.edit_b.Layout.Column=8;
btn = uibutton(gl, 'Text', 'Toggle TS sym', 'ButtonPushedFcn', @toggle_plot_style_callback);
btn.Layout.Row=2; btn.Layout.Column=9;

% ROW 3: Energy, Linewidth, Pressure
lbl = uilabel(gl, 'Text','RS E(mJ):'); lbl.Layout.Row=3; lbl.Layout.Column=1;
handles.edit_energy_r = uieditfield(gl, 'numeric', 'Value', params.energy_J_raman * 1000, 'ValueDisplayFormat', '%.2f'); handles.edit_energy_r.Layout.Row=3; handles.edit_energy_r.Layout.Column=2;
lbl = uilabel(gl, 'Text','TS E(mJ):'); lbl.Layout.Row=3; lbl.Layout.Column=3;
handles.edit_energy_t = uieditfield(gl, 'numeric', 'Value', params.energy_J_thomson * 1000, 'ValueDisplayFormat', '%.2f'); handles.edit_energy_t.Layout.Row=3; handles.edit_energy_t.Layout.Column=4;
lbl = uilabel(gl, 'Text','Linewidth:'); lbl.Layout.Row=3; lbl.Layout.Column=5;
handles.edit_linewidth = uieditfield(gl, 'numeric', 'Value', params.linewidth_initial, 'ValueDisplayFormat', '%.4f'); handles.edit_linewidth.Layout.Row=3; handles.edit_linewidth.Layout.Column=6;
lbl = uilabel(gl, 'Text','P(Pa):'); lbl.Layout.Row=3; lbl.Layout.Column=7;
handles.edit_pressure = uieditfield(gl, 'numeric', 'Value', params.pressure_Pa); handles.edit_pressure.Layout.Row=3; handles.edit_pressure.Layout.Column=[8,9];

% ROW 4: Temperature, Gas fractions, Recalibrate Raman button
lbl = uilabel(gl, 'Text','T(K):'); lbl.Layout.Row=4; lbl.Layout.Column=1;
handles.edit_temp_K = uieditfield(gl, 'numeric', 'Value', params.temp_K_initial,'Limits', [200, 400], 'ValueDisplayFormat', '%.2f'); handles.edit_temp_K.Layout.Row=4; handles.edit_temp_K.Layout.Column=2;
lbl = uilabel(gl, 'Text','N2:'); lbl.Layout.Row=4; lbl.Layout.Column=3;
handles.edit_n2_frac = uieditfield(gl, 'numeric', 'Value', params.mixing_ratios(1), 'Limits', [0, 1], 'ValueDisplayFormat', '%.3f'); handles.edit_n2_frac.Layout.Row=4; handles.edit_n2_frac.Layout.Column=4;
lbl = uilabel(gl, 'Text','O2:'); lbl.Layout.Row=4; lbl.Layout.Column=5;
handles.edit_o2_frac = uieditfield(gl, 'numeric', 'Value', params.mixing_ratios(2), 'Limits', [0, 1], 'ValueDisplayFormat', '%.3f'); handles.edit_o2_frac.Layout.Row=4; handles.edit_o2_frac.Layout.Column=6;
btn = uibutton(gl, 'Text', 'Recalib disp. from Raman', 'ButtonPushedFcn', @recalculate_calib_callback);
btn.Layout.Row=4; btn.Layout.Column=[7,9];

% ROW 5: Filters, Recalc Raman Spectrum button
handles.chk_remove_outliers = uicheckbox(gl, 'Text', 'RmOut', 'Value', false); 
handles.chk_remove_outliers.Layout.Row=5; handles.chk_remove_outliers.Layout.Column=[1,2];
handles.chk_median_filter = uicheckbox(gl, 'Text', 'MedF', 'Value', false); 
handles.chk_median_filter.Layout.Row=5; handles.chk_median_filter.Layout.Column=[3,4];
lbl = uilabel(gl, 'Text','Window:'); lbl.Layout.Row=5; lbl.Layout.Column=5;
handles.edit_medfilt_window = uieditfield(gl, 'numeric', 'Value', 5, 'Limits', [3, 15], 'RoundFractionalValues', 'on'); 
handles.edit_medfilt_window.Layout.Row=5; handles.edit_medfilt_window.Layout.Column=6;
btn = uibutton(gl, 'Text', 'Recalc Raman Spectrum', 'FontWeight', 'bold', 'BackgroundColor', [1, 0.9, 0.7], 'ButtonPushedFcn', @recalc_raman_spectrum_callback);
btn.Layout.Row=5; btn.Layout.Column=[7,9];

% ROW 6: Te, ne, alpha, gamma
lbl = uilabel(gl, 'Text','Te(eV):'); lbl.Layout.Row=6; lbl.Layout.Column=1;
handles.edit_Te_manual = uieditfield(gl, 'numeric', 'Value', data.deconv_results.Te_eV_corrected, 'Limits', [0.01, 50]); handles.edit_Te_manual.Layout.Row=6; handles.edit_Te_manual.Layout.Column=2;
lbl = uilabel(gl, 'Text','ne:'); lbl.Layout.Row=6; lbl.Layout.Column=3;
handles.edit_ne_manual = uieditfield(gl, 'numeric', 'Value', data.deconv_results.ne_corrected, 'ValueDisplayFormat','%.2e'); handles.edit_ne_manual.Layout.Row=6; handles.edit_ne_manual.Layout.Column=4;
lbl = uilabel(gl, 'Text','α:'); lbl.Layout.Row=6; lbl.Layout.Column=5;
alpha_initial = calculate_alpha_from_ne(data.deconv_results.Te_eV_corrected, data.deconv_results.ne_corrected, params.laser_wavelength);
handles.edit_alpha = uieditfield(gl, 'numeric', 'Value', alpha_initial, 'Limits', [1e-10, 5], 'Editable', 'off'); handles.edit_alpha.Layout.Row=6; handles.edit_alpha.Layout.Column=6;
lbl = uilabel(gl, 'Text','γ:'); lbl.Layout.Row=6; lbl.Layout.Column=7;


% ROW 7: Multi-starts and all buttons
handles.edit_gamma = uieditfield(gl, 'numeric', 'Value', 1.0, 'Limits', [0, 10]); handles.edit_gamma.Layout.Row=6; handles.edit_gamma.Layout.Column=8;
lbl = uilabel(gl, 'Text','Multi:'); lbl.Layout.Row=7; lbl.Layout.Column=1;
handles.edit_multistarts = uieditfield(gl, 'numeric', 'Value', 6, 'Limits', [1, 20], 'RoundFractionalValues', 'on'); 
handles.edit_multistarts.Layout.Row=7; handles.edit_multistarts.Layout.Column=2;
btn = uibutton(gl, 'Text', 'Recalc TS','BackgroundColor', [0.6, 1, 0.7], 'ButtonPushedFcn', @recalculate_fit_callback); 
btn.Layout.Row=7; btn.Layout.Column=3;
btn = uibutton(gl, 'Text', 'Shape Calib.', 'BackgroundColor', [0.3, 0.5, 0.9],'ButtonPushedFcn', @optimize_shape_callback); 
btn.Layout.Row=7; btn.Layout.Column=4;
btn = uibutton(gl, 'Text', 'Int.Calib.Gauss', 'FontWeight', 'bold', 'BackgroundColor', [0.9, 1, 0.9], 'ButtonPushedFcn', @reanalyze_all_callback); 
btn.Layout.Row=7; btn.Layout.Column=5;
btn = uibutton(gl, 'Text', 'Int.Calib.Full', 'FontWeight', 'bold', 'BackgroundColor', [0.85, 0.95, 1], 'ButtonPushedFcn', @intensity_calibration_callback); 
btn.Layout.Row=7; btn.Layout.Column=6;
btn = uibutton(gl, 'Text', '<', 'ButtonPushedFcn', @shift_fit_left_callback); 
btn.Layout.Row=7; btn.Layout.Column=8;
btn = uibutton(gl, 'Text', '>', 'ButtonPushedFcn', @shift_fit_right_callback); 
btn.Layout.Row=7; btn.Layout.Column=9;
handles.label_shift = uilabel(gl, 'Text', '0 px'); 
handles.label_shift.Layout.Row=7; handles.label_shift.Layout.Column=10;

% --- Utilities Panel - FINAL CLEAN LAYOUT (NO DUPLICATES) ---
util_gl = uigridlayout(utilities_panel);
util_gl.RowHeight = {'fit', 'fit', 'fit', 'fit', 'fit', 'fit'};
util_gl.ColumnWidth = {'fit', 'fit', 'fit', '1x', 'fit', 'fit'};

% ROW 1: Three utility buttons only
btn = uibutton(util_gl, 'Text', 'Decouple', 'ButtonPushedFcn', @decouple_plots_callback);
btn.Layout.Row = 1; btn.Layout.Column = 1;
btn = uibutton(util_gl, 'Text', 'EEDF', 'ButtonPushedFcn', @plot_eedf_callback);
btn.Layout.Row = 1; btn.Layout.Column = 2;
btn = uibutton(util_gl, 'Text', 'Save CSV', 'ButtonPushedFcn', @save_ts_data_callback);
btn.Layout.Row = 1; btn.Layout.Column = 3;

% ROW 2: Raman file
lbl = uilabel(util_gl, 'Text', 'Raman:', 'FontWeight', 'bold'); 
lbl.Layout.Row=2; lbl.Layout.Column=1;
handles.raman_path_display = uilabel(util_gl, 'Text', params.raman_file, 'Interpreter', 'none');
handles.raman_path_display.Layout.Row=2; handles.raman_path_display.Layout.Column=[2,4];
btn = uibutton(util_gl, 'Text', 'Load', 'ButtonPushedFcn', @load_raman_callback); 
btn.Layout.Row=2; btn.Layout.Column=5;
btn = uibutton(util_gl, 'Text', 'BG', 'ButtonPushedFcn', @load_raman_bg_callback); 
btn.Layout.Row=2; btn.Layout.Column=6;

% ROW 3: Raman BG
lbl = uilabel(util_gl, 'Text', 'Raman BG:', 'FontWeight', 'bold'); 
lbl.Layout.Row=3; lbl.Layout.Column=1;
handles.raman_bg_path_display = uilabel(util_gl, 'Text', 'No file loaded', 'Interpreter', 'none');
handles.raman_bg_path_display.Layout.Row=3; handles.raman_bg_path_display.Layout.Column=[2,6];

% ROW 4: Thomson file
lbl = uilabel(util_gl, 'Text', 'Thomson:', 'FontWeight', 'bold'); 
lbl.Layout.Row=4; lbl.Layout.Column=1;
handles.thomson_path_display = uilabel(util_gl, 'Text', params.thomson_file, 'Interpreter', 'none');
handles.thomson_path_display.Layout.Row=4; handles.thomson_path_display.Layout.Column=[2,4];
btn = uibutton(util_gl, 'Text', 'Load', 'ButtonPushedFcn', @load_thomson_callback); 
btn.Layout.Row=4; btn.Layout.Column=5;
btn = uibutton(util_gl, 'Text', 'BG', 'ButtonPushedFcn', @load_thomson_bg_callback); 
btn.Layout.Row=4; btn.Layout.Column=6;

% ROW 5: Thomson BG
lbl = uilabel(util_gl, 'Text', 'Thomson BG:', 'FontWeight', 'bold'); 
lbl.Layout.Row=5; lbl.Layout.Column=1;
handles.thomson_bg_path_display = uilabel(util_gl, 'Text', 'No file loaded', 'Interpreter', 'none');
handles.thomson_bg_path_display.Layout.Row=5; handles.thomson_bg_path_display.Layout.Column=[2,6];

% ROW 7: Status spanning all columns
handles.status_label = uilabel(util_gl, 'Text', 'Status: Ready', 'FontColor', [0, 0.5, 0], 'FontWeight', 'bold');
handles.status_label.Layout.Row = 6; handles.status_label.Layout.Column = [1,6];

populate_file_list(handles);
guidata(fig, handles);
fig.SizeChangedFcn = @(src, event) handle_resize(src, handles);
update_images_and_status(handles);

    function populate_file_list(handles)
        [file_path, ~, ~] = fileparts(handles.thomson_file_original);

        if isempty(file_path) || ~isfolder(file_path)
            handles.file_listbox.Items = {'No valid folder'};
            handles.file_listbox.ItemsData = {''};
            return;
        end

        % Get all image files
        image_extensions = {'*.tif', '*.tiff', '*.png', '*.jpg', '*.jpeg'};
        all_files = [];

        for i = 1:length(image_extensions)
            files = dir(fullfile(file_path, image_extensions{i}));
            all_files = [all_files; files];
        end

        if isempty(all_files)
            handles.file_listbox.Items = {'No image files found'};
            handles.file_listbox.ItemsData = {''};
            return;
        end

        % Sort by date (most recent first)
        [~, idx] = sort([all_files.datenum], 'descend');
        all_files = all_files(idx);

        % Populate listbox
        file_names = {all_files.name}';
        file_paths = cellfun(@(x) fullfile(file_path, x), file_names, 'UniformOutput', false);

        handles.file_listbox.Items = file_names;
        handles.file_listbox.ItemsData = file_paths;

        % Highlight current file
        [~, current_name, current_ext] = fileparts(handles.thomson_file_original);
        current_full_name = [current_name, current_ext];
        current_idx = find(strcmp(file_names, current_full_name));
        if ~isempty(current_idx)
            handles.file_listbox.Value = file_paths{current_idx};
        end
    end

    function file_selected_callback(hObject, ~)
        handles = guidata(hObject);

        selected_file = handles.file_listbox.Value;

        if isempty(selected_file) || strcmp(selected_file, '')
            return;
        end

        if strcmp(selected_file, handles.thomson_file_original)
            return;
        end

        update_status(handles, 'Loading new Thomson file...', [0, 0, 1]);

        handles.thomson_file_original = selected_file;
        handles = update_images_and_status(handles);

        current_params = handles.params;
        current_params.thomson_file = selected_file;
        current_params.thomson_img_raw = double(imread(selected_file));
        if ~isempty(handles.thomson_bg_file)
            current_params.thomson_img_raw = current_params.thomson_img_raw - double(imread(handles.thomson_bg_file));
        end
        current_params.yc = handles.edit_yc.Value;
        current_params.wy = handles.edit_wy.Value;
        current_params.exclusion_cols = handles.edit_start.Value : handles.edit_end.Value;
        current_params.shots_thomson = handles.edit_shots_t.Value;
        current_params.energy_J_thomson = handles.edit_energy_t.Value / 1000;
        current_params.shots_raman = handles.edit_shots_r.Value;
        current_params.energy_J_raman = handles.edit_energy_r.Value / 1000;
        current_params.signal_rows = (current_params.yc - floor(current_params.wy/2)):(current_params.yc + floor(current_params.wy/2));
        current_params.temp_K_initial = handles.edit_temp_K.Value;
        current_params.mixing_ratios = [handles.edit_n2_frac.Value, handles.edit_o2_frac.Value];

        % Sync params with GUI before analysis
        [new_results, new_plot_data] = run_full_analysis(current_params);

        handles.params = current_params;
        handles.plot_data = new_plot_data;

        thomson_spectrum = apply_ts_filters(new_plot_data.thomson_spectrum_processed, ...
            handles.chk_remove_outliers.Value, ...
            handles.chk_median_filter.Value, ...
            handles.edit_medfilt_window.Value);

        handles.current_ts_spectrum = new_plot_data.thomson_spectrum_processed;
        handles.current_results.uncorrected = new_results.uncorrected;
        handles.current_deconv = new_results.deconv_results;

        handles.edit_Te_manual.Value = new_plot_data.deconv_results.Te_eV_corrected;
        handles.edit_ne_manual.Value = new_plot_data.deconv_results.ne_corrected;
        alpha_new = calculate_alpha_from_ne(new_plot_data.deconv_results.Te_eV_corrected, new_plot_data.deconv_results.ne_corrected, current_params.laser_wavelength);
        try
            handles.edit_alpha.Value = alpha_new;
        catch ME
            % If value is out of range, use a default or clamp it
            alpha_limits = handles.edit_alpha.Limits;
            if alpha_new < alpha_limits(1)
                handles.edit_alpha.Value = alpha_limits(1);
                warning('Alpha too small (%.4e), set to minimum limit %.4e', alpha_new, alpha_limits(1));
            elseif alpha_new > alpha_limits(2)
                handles.edit_alpha.Value = alpha_limits(2);
                warning('Alpha too large (%.4e), set to maximum limit %.4e', alpha_new, alpha_limits(2));
            else
                rethrow(ME); % Some other error
            end
        end

        handles.edit_gamma.Value = 1.0;

        update_plots(handles);

        fprintf('Switched to file: %s\n', selected_file);
        update_status(handles, 'Ready', [0, 0.5, 0]);
        guidata(hObject, handles);
    end

    function refresh_file_list_callback(hObject, ~)
        handles = guidata(hObject);
        populate_file_list(handles);
        guidata(hObject, handles);
    end

    function refresh_images_roi_callback(hObject, ~)
        handles = guidata(hObject);
        update_status(handles, 'Refreshing images and ROIs...', [0, 0, 1]);
        disp('--- Refreshing Images and ROI Overlays ---');

        yc = handles.edit_yc.Value;
        wy = handles.edit_wy.Value;
        excl_start = handles.edit_start.Value;
        excl_end = handles.edit_end.Value;

        raman_raw = handles.params.raman_img_raw;
        thomson_raw = handles.params.thomson_img_raw;
        [nRows, nCols] = size(raman_raw);

        params_adjusted = false;
        if (yc + floor(wy/2)) > nRows || (yc - floor(wy/2)) < 1
            yc = round(nRows / 2);
            wy = min(wy, floor(nRows/2));
            if mod(wy, 2) ~= 0, wy = wy - 1; end
            handles.edit_yc.Value = yc;
            handles.edit_wy.Value = wy;
            params_adjusted = true;
        end

        if excl_start < 1 || excl_end > nCols || excl_start >= excl_end
            new_center_col = round(nCols / 2);
            new_excl_width = min(excl_end - excl_start + 1, round(nCols / 10));
            excl_start = new_center_col - floor(new_excl_width / 2);
            excl_end = new_center_col + floor(new_excl_width / 2);
            handles.edit_start.Value = excl_start;
            handles.edit_end.Value = excl_end;
            params_adjusted = true;
        end

        if params_adjusted
            fprintf('INFO: ROI/Exclusion parameters were adjusted to fit image dimensions.\n');
        end

        signal_rows = (yc - floor(wy/2)):(yc + floor(wy/2));
        exclusion_cols = excl_start:excl_end;

        handles.params.yc = yc;
        handles.params.wy = wy;
        handles.params.signal_rows = signal_rows;
        handles.params.exclusion_cols = exclusion_cols;

        raman_processed = block_central_pixels(raman_raw, signal_rows, exclusion_cols);
        thomson_processed = block_central_pixels(thomson_raw, signal_rows, exclusion_cols);

        set(handles.h_raman_img, 'CData', medfilt2(raman_processed, [8 8]));
        set(handles.h_thomson_img, 'CData', medfilt2(thomson_processed, [8 8]));

        set(handles.h_raman_roi_rect, 'Position', [1, signal_rows(1), nCols-1, length(signal_rows)]);
        set(handles.h_raman_excl_rect, 'Position', [excl_start, 1, length(exclusion_cols), nRows]);
        set(handles.h_thomson_roi_rect, 'Position', [1, signal_rows(1), nCols-1, length(signal_rows)]);
        set(handles.h_excl_rect, 'Position', [excl_start, 1, length(exclusion_cols), nRows]);

        pixels_exp = 1:nCols;
        wavelength_axis = handles.edit_a.Value * pixels_exp + handles.edit_b.Value;

        experimental_raman = mean(raman_processed(signal_rows, :), 1);
        experimental_raman = subtract_wing_baseline(experimental_raman);
        handles.plot_data.raman_results.experimental_spectrum = experimental_raman;
        handles.plot_data.raman_results.wavelength_axis = wavelength_axis;

        set(handles.h_raman_spec, 'XData', wavelength_axis, 'YData', experimental_raman);

        if isfield(handles.plot_data.raman_results, 'theoretical_spectrum_model')
            model_scaled = handles.plot_data.raman_results.theoretical_spectrum_model * ...
                (max(experimental_raman) / max(handles.plot_data.raman_results.theoretical_spectrum_model));
            set(handles.h_raman_model, 'XData', wavelength_axis, 'YData', model_scaled);
        end

        xlim(handles.ax3, [min(wavelength_axis), max(wavelength_axis)]);

        raw_thomson = mean(thomson_processed(signal_rows, :), 1);
        thomson_spectrum = subtract_wing_baseline(raw_thomson);

        thomson_spectrum = apply_ts_filters(thomson_spectrum, ...
            handles.chk_remove_outliers.Value, ...
            handles.chk_median_filter.Value, ...
            handles.edit_medfilt_window.Value);

        handles.current_ts_spectrum = thomson_spectrum;
        handles.plot_data.thomson_spectrum_processed = thomson_spectrum;
        handles.plot_data.wavelength_nm = wavelength_axis;

        fit_mask = ~ismember(1:length(wavelength_axis), exclusion_cols);
        if isfield(handles, 'h_data_points') && isvalid(handles.h_data_points)
            set(handles.h_data_points, 'XData', wavelength_axis(fit_mask), 'YData', thomson_spectrum(fit_mask));
        end

        handles.plot_data.raman_img_processed = raman_processed;
        handles.plot_data.thomson_img_processed = thomson_processed;

        axis(handles.ax1, 'image');
        axis(handles.ax2, 'image');

        fprintf('Images and ROIs refreshed: yc=%d, wy=%d, excl=[%d:%d]\n', yc, wy, excl_start, excl_end);
        update_status(handles, 'Ready', [0, 0.5, 0]);
        guidata(hObject, handles);
    end

    function recalc_raman_spectrum_callback(hObject, ~)
        handles = guidata(hObject);
        update_status(handles, 'Recalculating Raman spectrum...', [0, 0, 1]);
        disp('--- Recalculating Raman Spectrum from GUI Parameters ---');

        % Get current values from GUI
        a_coeff = handles.edit_a.Value;
        b_coeff = handles.edit_b.Value;
        linewidth = handles.edit_linewidth.Value;
        pressure_Pa = handles.edit_pressure.Value;
        temp_K = handles.edit_temp_K.Value;
        n2_frac = handles.edit_n2_frac.Value;
        o2_frac = handles.edit_o2_frac.Value;

        % Validate that fractions sum to 1 (or close to it)
        total_frac = n2_frac + o2_frac;
        if abs(total_frac - 1.0) > 0.01
            uialert(ancestor(hObject, 'figure'), ...
                sprintf('N2 and O2 fractions must sum to 1.0 (current sum: %.3f)', total_frac), ...
                'Invalid Gas Composition');
            update_status(handles, 'Error: Invalid gas fractions', [1, 0, 0]);
            return;
        end

        % Normalize if close but not exact
        if abs(total_frac - 1.0) > 1e-6
            n2_frac = n2_frac / total_frac;
            o2_frac = o2_frac / total_frac;
            handles.edit_n2_frac.Value = n2_frac;
            handles.edit_o2_frac.Value = o2_frac;
        end

        % Calculate new wavelength axis
        pixels_exp = 1:size(handles.params.raman_img_raw, 2);
        wavelength_axis = a_coeff * pixels_exp + b_coeff;

        % Calculate new theoretical spectrum with updated parameters
        species_list = {'N2', 'O2'};
        mixing_ratios = [n2_frac, o2_frac];

        theoretical_spectrum = calculate_theoretical_spectrum(wavelength_axis, ...
            species_list, mixing_ratios, pressure_Pa, temp_K, linewidth, ...
            handles.params.laser_wavelength);

        % Recalculate calibration factor
        experimental_integral = trapz(wavelength_axis, ...
            handles.plot_data.raman_results.experimental_spectrum);
        theoretical_integral = trapz(wavelength_axis, theoretical_spectrum);
        total_energy_raman = (handles.edit_energy_r.Value/1000) * handles.edit_shots_r.Value;

        if theoretical_integral > 0
            new_calibration_factor = experimental_integral / ...
                (total_energy_raman * theoretical_integral);
        else
            new_calibration_factor = 1e-15;
        end

        % Update stored data
        handles.plot_data.raman_results.wavelength_axis = wavelength_axis;
        handles.plot_data.raman_results.theoretical_spectrum_model = theoretical_spectrum;
        handles.plot_data.raman_results.calibration_factor = new_calibration_factor;
        handles.temp_K_current = temp_K;
        handles.n2_fraction = n2_frac;
        handles.o2_fraction = o2_frac;

        % Update Raman spectrum plot
        set(handles.h_raman_spec, 'XData', wavelength_axis, ...
            'YData', handles.plot_data.raman_results.experimental_spectrum);

        model_scaled = theoretical_spectrum * ...
            (max(handles.plot_data.raman_results.experimental_spectrum) / ...
            max(theoretical_spectrum));
        set(handles.h_raman_model, 'XData', wavelength_axis, 'YData', model_scaled);

        set(handles.h_calib_title, 'String', ...
            sprintf('Raman Spectrum (Calib. Factor: %.2e)\nT=%.1f K, N2=%.1f%%, O2=%.1f%%', ...
            new_calibration_factor, temp_K, n2_frac*100, o2_frac*100));

        xlim(handles.ax3, [min(wavelength_axis), max(wavelength_axis)]);

        fprintf('Raman recalculated: Calib=%.2e, T=%.1f K, N2=%.3f, O2=%.3f, LW=%.4f nm\n', ...
            new_calibration_factor, temp_K, n2_frac, o2_frac, linewidth);

        update_status(handles, 'Ready', [0, 0.5, 0]);
        guidata(hObject, handles);
    end


    function save_axis_plot(ax, plot_type)
        % Get the current handles structure to ensure up-to-date info
        fig = ancestor(ax, 'figure');
        handles = guidata(fig);

        % Determine which file to use as the base name
        if contains(plot_type, 'raman')
            reference_file = handles.raman_file_original;
        else % 'thomson' or 'ts'
            reference_file = handles.thomson_file_original;
        end

        [file_path, file_name, ~] = fileparts(reference_file);

        % Build filename based on plot type
        switch plot_type
            case 'raman_image'
                save_name = fullfile(file_path, sprintf('%s_RamanImage.jpg', file_name));
            case 'thomson_image'
                save_name = fullfile(file_path, sprintf('%s_ThomsonImage.jpg', file_name));
            case 'raman_spectrum'
                calib = handles.plot_data.raman_results.calibration_factor;
                save_name = fullfile(file_path, sprintf('%s_RamanSpectrum_CalibFactor_%.2e.jpg', file_name, calib));
            case 'ts_spectrum'
                Te = handles.edit_Te_manual.Value;
                ne = handles.edit_ne_manual.Value;

                % Determine which fit is currently visible
                fit_type = 'Unknown';
                ssr = NaN;
                if isfield(handles, 'h_coherent_fit') && isvalid(handles.h_coherent_fit) && strcmp(get(handles.h_coherent_fit, 'Visible'), 'on')
                    displayName = get(handles.h_coherent_fit, 'DisplayName');
                    if contains(displayName, 'Intensity Cal')
                        fit_type = 'IntensityCal';
                        if isfield(handles, 'intensity_fit_ssr'), ssr = handles.intensity_fit_ssr; end
                    elseif contains(displayName, 'Shape Fit')
                        fit_type = 'ShapeFit';
                        if isfield(handles, 'shape_fit_ssr'), ssr = handles.shape_fit_ssr; end
                        % MODIFIED: Added case for 'Manual Fit'
                    elseif contains(displayName, 'Manual Fit')
                        fit_type = 'ManualFit';
                        ssr = NaN; % Manual calculations don't have an SSR
                    else
                        fit_type = 'CoherentFit';
                    end
                else
                    fit_type = 'Gaussian';
                    ssr = handles.plot_data.ssr_gauss;
                end

                if isfinite(ssr)
                    save_name = fullfile(file_path, sprintf('%s_TSSpectrum_%s_Te%.2feV_ne%.2e_SSR%.2e.jpg', file_name, fit_type, Te, ne, ssr));
                else
                    save_name = fullfile(file_path, sprintf('%s_TSSpectrum_%s_Te%.2feV_ne%.2e.jpg', file_name, fit_type, Te, ne));
                end
            otherwise
                save_name = fullfile(file_path, sprintf('%s.jpg', file_name));
        end

        temp_fig = figure('Visible', 'off');
        temp_ax = copyobj(ax, temp_fig);
        set(temp_ax, 'Units', 'normalized', 'Position', [0.1, 0.1, 0.8, 0.8]);

        saveas(temp_fig, save_name);
        close(temp_fig);

        fprintf('Saved plot to: %s\n', save_name);
        uialert(fig, sprintf('Plot saved to:\n%s', save_name), 'Save Successful');
    end

    function handle_resize(fig, handles)
        % This function handles dynamic resizing of the GUI
        % No need to do anything - uigridlayout handles it automatically
        % This is just a placeholder for future custom resize behavior if needed
        drawnow;
    end

    function update_status(handles, message, color)
        if nargin < 3, color = [0, 0, 0]; end
        handles.status_label.Text = ['Status: ' message];
        handles.status_label.FontColor = color;
        drawnow;
    end

    function load_raman_callback(hObject, ~)
        handles = guidata(hObject);
        [start_path, ~, ~] = fileparts(handles.raman_file_original);
        [file, path] = uigetfile({'*.tif;*.tiff';'*.png';'*.*'}, 'Select Raman Image', start_path);
        if file ~= 0
            handles.raman_file_original = fullfile(path, file);
            handles = update_images_and_status(handles);
            fprintf('Loaded Raman image: %s\n', handles.raman_file_original);
        end
        guidata(hObject, handles);
    end

    function load_raman_bg_callback(hObject, ~)
        handles = guidata(hObject);
        [start_path, ~, ~] = fileparts(handles.raman_file_original);
        [file, path] = uigetfile({'*.tif;*.tiff';'*.png';'*.*'}, 'Select Raman Background', start_path);
        if file ~= 0
            handles.raman_bg_file = fullfile(path, file);
            handles = update_images_and_status(handles);
            fprintf('Loaded Raman background: %s\n', handles.raman_bg_file);
        end
        guidata(hObject, handles);
    end

    function load_thomson_callback(hObject, ~)
        handles = guidata(hObject);
        [start_path, ~, ~] = fileparts(handles.thomson_file_original);
        [file, path] = uigetfile({'*.tif;*.tiff';'*.png';'*.*'}, 'Select Thomson Image', start_path);
        if file ~= 0
            handles.thomson_file_original = fullfile(path, file);
            handles = update_images_and_status(handles);
            fprintf('Loaded Thomson image: %s\n', handles.thomson_file_original);
        end
        guidata(hObject, handles);
    end

    function load_thomson_bg_callback(hObject, ~)
        handles = guidata(hObject);
        [start_path, ~, ~] = fileparts(handles.thomson_file_original);
        [file, path] = uigetfile({'*.tif;*.tiff';'*.png';'*.*'}, 'Select Thomson Background', start_path);
        if file ~= 0
            handles.thomson_bg_file = fullfile(path, file);
            handles = update_images_and_status(handles);
            fprintf('Loaded Thomson background: %s\n', handles.thomson_bg_file);
        end
        guidata(hObject, handles);
    end

    function handles = update_images_and_status(handles)
        raman_raw = double(imread(handles.raman_file_original));
        if ~isempty(handles.raman_bg_file), raman_raw = raman_raw - double(imread(handles.raman_bg_file)); end
        thomson_raw = double(imread(handles.thomson_file_original));
        if ~isempty(handles.thomson_bg_file), thomson_raw = thomson_raw - double(imread(handles.thomson_bg_file)); end

        [nRows, nCols] = size(raman_raw);
        params_adjusted = false;

        yc = handles.edit_yc.Value;
        wy = handles.edit_wy.Value;
        if (yc + floor(wy/2)) > nRows || (yc - floor(wy/2)) < 1
            handles.edit_yc.Value = round(nRows / 2);
            handles.edit_wy.Value = min(wy, floor(nRows/2));
            params_adjusted = true;
        end

        excl_start = handles.edit_start.Value;
        excl_end = handles.edit_end.Value;
        if excl_start < 1 || excl_end > nCols
            new_center_col = round(nCols / 2);
            new_excl_width = min(excl_end - excl_start + 1, round(nCols / 10));
            handles.edit_start.Value = new_center_col - floor(new_excl_width / 2);
            handles.edit_end.Value = new_center_col + floor(new_excl_width / 2);
            params_adjusted = true;
        end

        if params_adjusted, fprintf('INFO: ROI and/or Exclusion parameters were auto-adjusted to fit new image dimensions.\n'); end
        drawnow;

        handles.params.raman_img_raw = raman_raw;
        handles.params.thomson_img_raw = thomson_raw;

        current_params = handles.params;
        current_params.signal_rows = (handles.edit_yc.Value - floor(handles.edit_wy.Value/2)):(handles.edit_yc.Value + floor(handles.edit_wy.Value/2));
        current_params.exclusion_cols = handles.edit_start.Value : handles.edit_end.Value;
        raman_processed = block_central_pixels(raman_raw, current_params.signal_rows, current_params.exclusion_cols);
        thomson_processed = block_central_pixels(thomson_raw, current_params.signal_rows, current_params.exclusion_cols);
        set(handles.h_raman_img, 'CData', medfilt2(raman_processed,[8 8]));
        set(handles.h_thomson_img, 'CData', medfilt2(thomson_processed, [8 8]));

        axis(handles.ax1, 'image');% Raman now
        axis(handles.ax2, 'image');% Thomson now

        % Corrected title logic for ax1 (Raman)
        if ~isempty(handles.raman_bg_file)
            title(handles.ax1, 'Raman Image (BG Subtracted)');
        else
            title(handles.ax1, 'Raman Image (Raw)');
        end

        % Corrected title logic for ax2 (Thomson)
        if ~isempty(handles.thomson_bg_file)
            title(handles.ax2, 'Thomson Image (BG Subtracted)');
        else
            title(handles.ax2, 'Thomson Image (Raw)');
        end

        handles.raman_path_display.Text = handles.raman_file_original;
        handles.thomson_path_display.Text = handles.thomson_file_original;
        if ~isempty(handles.raman_bg_file), handles.raman_bg_path_display.Text = handles.raman_bg_file; handles.raman_bg_path_display.FontColor = [0 0.4470 0.7410];
        else, handles.raman_bg_path_display.Text = 'No file loaded'; handles.raman_bg_path_display.FontColor = [0.6 0.6 0.6]; end
        if ~isempty(handles.thomson_bg_file), handles.thomson_bg_path_display.Text = handles.thomson_bg_file; handles.thomson_bg_path_display.FontColor = [0 0.4470 0.7410];
        else, handles.thomson_bg_path_display.Text = 'No file loaded'; handles.thomson_bg_path_display.FontColor = [0.6 0.6 0.6]; end
    end

    function recalculate_fit_callback(hObject, ~)
        handles = guidata(hObject);
        update_status(handles, 'Recalculating TS spectrum...', [0, 0, 1]);
        disp('--- Recalculating TS spectrum from manual Te, ne, alpha, and gamma ---');

        % Read values from GUI
        manual_Te = handles.edit_Te_manual.Value;
        manual_ne = handles.edit_ne_manual.Value;
        manual_alpha = handles.edit_alpha.Value;
        manual_gamma = handles.edit_gamma.Value;

        % Recalculate alpha from Te and ne to keep consistent
        laser_wl = handles.params.laser_wavelength;
        manual_alpha = calculate_alpha_from_ne(manual_Te, manual_ne, laser_wl);
        handles.edit_alpha.Value = manual_alpha;  % Update the alpha field


        % Get necessary parameters
        laser_wl = handles.params.laser_wavelength;
        wavelength_nm = handles.plot_data.wavelength_nm;
        calib_factor = handles.plot_data.raman_results.calibration_factor;
        total_energy_thomson = (handles.edit_energy_t.Value/1000) * handles.edit_shots_t.Value;
        sigma_T_diff = 7.94e-30;

        % Generate the coherent theoretical spectrum shape
        shape = calculate_coherent_ts(wavelength_nm, manual_Te, manual_alpha, laser_wl);
        sum_shape = sum(shape);
        if sum_shape > 1e-12
            normalized_shape = shape / sum_shape;
        else
            normalized_shape = zeros(size(wavelength_nm));
        end

        % Apply gamma as duty cycle and calculate total signal
        total_signal = calib_factor * total_energy_thomson * sigma_T_diff * manual_ne * manual_gamma;
        updated_spectrum = total_signal * normalized_shape;

        % Update or create the plot
        if isfield(handles, 'h_coherent_fit') && isvalid(handles.h_coherent_fit)
            % MODIFIED: Update DisplayName to 'Manual Fit'
            set(handles.h_coherent_fit, 'YData', updated_spectrum, 'Visible', 'on', 'DisplayName', 'Manual Fit');
            if isfield(handles, 'h_gauss_fit') && isvalid(handles.h_gauss_fit)
                set(handles.h_gauss_fit, 'Visible', 'off');
            end
        else
            cla(handles.ax4);
            fit_mask = ~ismember(1:length(wavelength_nm), handles.edit_start.Value:handles.edit_end.Value);
            handles.h_data_points = plot(handles.ax4, wavelength_nm(fit_mask), handles.current_ts_spectrum(fit_mask), 'b-', 'DisplayName', 'Data');
            hold(handles.ax4, 'on');
            % MODIFIED: Set DisplayName to 'Manual Fit'
            handles.h_coherent_fit = plot(handles.ax4, wavelength_nm, updated_spectrum, 'g-', 'LineWidth', 2, 'DisplayName', 'Manual Fit');
            grid(handles.ax4, 'on'); legend(handles.ax4, 'Location', 'best');
            xlabel(handles.ax4, 'Wavelength (nm)'); ylabel(handles.ax4, 'Intensity (counts)');
        end

        % MODIFIED: Update title to 'Manual Fit'
        title(handles.ax4, sprintf('Manual Fit: Te=%.2f eV, n_{e}=%.2e m^{-3}, \\alpha=%.2f, \\gamma=%.2f', manual_Te, manual_ne, manual_alpha, manual_gamma));

        update_status(handles, 'Ready', [0, 0.5, 0]);
        guidata(hObject, handles);
    end

    function intensity_calibration_callback(hObject, ~)
        handles = guidata(hObject);
        update_status(handles, 'Performing intensity calibration (2-stage fit)...', [1, 0.5, 0]);
        disp('--- Intensity Calibration: Stage 1 - Shape fitting ---');

        wavelength_nm = handles.plot_data.wavelength_nm;
        thomson_spectrum = handles.current_ts_spectrum;
        exclusion_cols = handles.edit_start.Value:handles.edit_end.Value;
        fit_mask = ~ismember(1:length(wavelength_nm), exclusion_cols);
        laser_wl = handles.params.laser_wavelength;

        % === NEW: Create wavelength-dependent weights centered at 532 nm ===
        central_wl = 532;  % nm - prioritize this region
        weight_sigma = 5;  % nm - controls how quickly weights decay from center
        weights_full = exp(-((wavelength_nm - central_wl).^2) / (2 * weight_sigma^2));
        weights_full = weights_full / max(weights_full);  % Normalize to [0, 1]
        weights = weights_full(fit_mask);  % Extract only the fit region weights
        fprintf('Using Gaussian weights centered at %.1f nm (sigma=%.1f nm) to prioritize central region\n', central_wl, weight_sigma);

        shape_objective = @(p, x) coherent_fitter_direct(p, x, laser_wl);

        Te_initial = handles.edit_Te_manual.Value;
        ne_initial = handles.edit_ne_manual.Value;
        alpha_initial = calculate_alpha_from_ne(Te_initial, ne_initial, laser_wl);
        amp_initial = max(thomson_spectrum(fit_mask));

        fprintf('Starting intensity calibration from: Te=%.2f eV, ne=%.2e m-3, alpha=%.2f\n', Te_initial, ne_initial, alpha_initial);

        %%%%%%%%%% MultiStart fitting
        num_starts = handles.edit_multistarts.Value;
        best_p_fit = [];
        best_resnorm = inf;
        best_jacobian = [];
        options = optimoptions('lsqcurvefit', 'Display', 'final'); % Display set to 'final'

        for start_idx = 1:num_starts
            if start_idx == 1
                p0_current = [Te_initial, alpha_initial, amp_initial];
            else
                p0_current = [Te_initial*(0.8 + 0.4*rand()), alpha_initial*(0.5 + 1.0*rand()), amp_initial*(0.8 + 0.4*rand())];
            end

            lb = [0.1, 0.01, 0];
            ub = [15.0, 5.0, 3 * amp_initial];

            try
                % === NEW: Use weighted fitting - multiply both data and model by sqrt(weights) ===
                sqrt_weights = sqrt(weights);
                weighted_data = thomson_spectrum(fit_mask) .* sqrt_weights;
                weighted_objective = @(p, x) shape_objective(p, x) .* sqrt_weights;

                [p_fit_current, resnorm_current, ~, exitflag, ~, ~, jacobian_current] = lsqcurvefit(weighted_objective, p0_current, wavelength_nm(fit_mask), weighted_data, lb, ub, options);

                if exitflag > 0 && resnorm_current < best_resnorm
                    best_p_fit = p_fit_current;
                    best_resnorm = resnorm_current;
                    best_jacobian = jacobian_current;
                    fprintf('  Intensity cal multi-start %d: Te=%.2f, alpha=%.2f, SSR=%.2e\n', start_idx, p_fit_current(1), p_fit_current(2), resnorm_current);
                end
            catch ME
                fprintf('  Multi-start %d failed: %s\n', start_idx, ME.message);
            end
        end

        if isempty(best_p_fit)
            uialert(ancestor(hObject, 'figure'), 'Intensity calibration fit failed to converge.', 'Error');
            update_status(handles, 'Fit failed', [1, 0, 0]);
            return;
        end

        p_fit = best_p_fit;

        Te_fit = p_fit(1);
        alpha_fit = p_fit(2);
        amp_fit = p_fit(3);

        disp('--- Intensity Calibration: Stage 2 - Calculate ne and realistic errors via Monte Carlo ---');

        fitted_shape_spectrum = coherent_fitter_direct(p_fit, wavelength_nm, laser_wl);
        total_counts = sum(fitted_shape_spectrum);

        calib_factor = handles.plot_data.raman_results.calibration_factor;
        total_energy_thomson = (handles.edit_energy_t.Value/1000) * handles.edit_shots_t.Value;
        sigma_T_diff = 7.94e-30;

        ne_fit = total_counts / (calib_factor * total_energy_thomson * sigma_T_diff);

        % === NEW: Monte Carlo error estimation - find parameter ranges that produce visually acceptable fits ===
        fprintf('Calculating realistic parameter errors via Monte Carlo sampling...\n');

        % Define search ranges (±30% to capture typical variations)
        Te_range = linspace(Te_fit * 0.7, Te_fit * 1.3, 25);
        ne_range = linspace(ne_fit * 0.7, ne_fit * 1.3, 25);

        % Define acceptable SSR threshold (1.5x best fit = visually indistinguishable)
        % Adjust this multiplier (1.5) to be more/less conservative
        acceptable_ssr_threshold = best_resnorm * 1.5;
        fprintf('  SSR threshold for acceptable fits: %.2e (best fit SSR = %.2e)\n', acceptable_ssr_threshold, best_resnorm);

        acceptable_Te = [];
        acceptable_ne = [];

        % Grid search over parameter space
        for Te_test = Te_range
            for ne_test = ne_range
                % Calculate alpha for this Te/ne combination
                alpha_test = calculate_alpha_from_ne(Te_test, ne_test, laser_wl);

                % Generate synthetic spectrum with these parameters
                test_spectrum = coherent_fitter_direct([Te_test, alpha_test, amp_fit], wavelength_nm, laser_wl);

                % Calculate weighted residuals (same weights as fitting)
                sqrt_weights_full = sqrt(weights_full(fit_mask));
                test_residuals_weighted = (thomson_spectrum(fit_mask) - test_spectrum(fit_mask)) .* sqrt_weights_full;
                test_ssr = sum(test_residuals_weighted.^2);

                % If this fit is within acceptable threshold, add to list
                if test_ssr < acceptable_ssr_threshold
                    acceptable_Te = [acceptable_Te, Te_test];
                    acceptable_ne = [acceptable_ne, ne_test];
                end
            end
        end

        % Calculate errors as half the range of acceptable values
        if ~isempty(acceptable_Te)
            Te_fit_err = (max(acceptable_Te) - min(acceptable_Te)) / 2;
            ne_fit_err = (max(acceptable_ne) - min(acceptable_ne)) / 2;
            fprintf('  Monte Carlo results: %d parameter sets within acceptable SSR\n', length(acceptable_Te));
            fprintf('  Te range: %.2f to %.2f eV (error = %.2f eV, %.1f%%)\n', min(acceptable_Te), max(acceptable_Te), Te_fit_err, 100*Te_fit_err/Te_fit);
            fprintf('  ne range: %.2e to %.2e m^-3 (error = %.2e m^-3, %.1f%%)\n', min(acceptable_ne), max(acceptable_ne), ne_fit_err, 100*ne_fit_err/ne_fit);
        else
            % Fallback: use 15% as empirically observed typical error
            Te_fit_err = Te_fit * 0.15;
            ne_fit_err = ne_fit * 0.15;
            fprintf('  Warning: No parameter sets found within threshold. Using fallback 15%% errors.\n');
        end

        residuals = thomson_spectrum(fit_mask) - fitted_shape_spectrum(fit_mask);
        ssr_intensity = sum(residuals.^2);
        tss = sum((thomson_spectrum(fit_mask) - mean(thomson_spectrum(fit_mask))).^2);
        r_squared_intensity = 1 - ssr_intensity/tss;

        fprintf('Intensity Cal Results: Te=%.2f +/- %.2f eV, alpha=%.2f, ne=%.2e +/- %.2e m-3, SSR=%.2e, R²=%.3f\n', Te_fit, Te_fit_err, alpha_fit, ne_fit, ne_fit_err, ssr_intensity, r_squared_intensity);

        handles.intensity_fit_ssr = ssr_intensity;
        handles.intensity_fit_r2 = r_squared_intensity;

        handles.edit_Te_manual.Value = Te_fit;
        handles.edit_alpha.Value = alpha_fit;
        handles.edit_ne_manual.Value = ne_fit;
        handles.edit_gamma.Value = 1.0;

        if isfield(handles, 'h_coherent_fit') && isvalid(handles.h_coherent_fit)
            set(handles.h_coherent_fit, 'YData', fitted_shape_spectrum, 'DisplayName', 'Intensity Cal', 'Visible', 'on');
            if isfield(handles, 'h_gauss_fit') && isvalid(handles.h_gauss_fit), set(handles.h_gauss_fit, 'Visible','off'); end
        else
            cla(handles.ax4);
            handles.h_data_points = plot(handles.ax4, wavelength_nm(fit_mask), thomson_spectrum(fit_mask), 'b-', 'DisplayName', 'Data');
            hold(handles.ax4, 'on');
            handles.h_coherent_fit = plot(handles.ax4, wavelength_nm, fitted_shape_spectrum, 'g-', 'LineWidth', 2, 'DisplayName', 'Intensity Cal');
            grid(handles.ax4, 'on'); legend(handles.ax4, 'Location', 'best');
            xlabel(handles.ax4, 'Wavelength (nm)'); ylabel(handles.ax4, 'Intensity (counts)');
        end

        title(handles.ax4, sprintf('Intensity Cal (weighted fit): Te=%.2f \\pm %.2f eV, n_{e}=%.2e \\pm %.2e m^{-3}, \\alpha=%.2f\nSSR=%.2e, R^2=%.3f', Te_fit, Te_fit_err, ne_fit, ne_fit_err, alpha_fit, ssr_intensity, r_squared_intensity));

        update_status(handles, 'Ready', [0, 0.5, 0]);
        guidata(hObject, handles);
    end

    function optimize_shape_callback(hObject, ~)
        handles = guidata(hObject);
        update_status(handles, 'Optimizing shape fit...', [1, 0.5, 0]);
        disp('--- Performing multi-start coherent shape fit ---');

        wavelength_nm = handles.plot_data.wavelength_nm;
        thomson_spectrum = handles.current_ts_spectrum;
        exclusion_cols = handles.edit_start.Value:handles.edit_end.Value;
        fit_mask = ~ismember(1:length(wavelength_nm), exclusion_cols);
        laser_wl = handles.params.laser_wavelength;

        Te_initial = handles.edit_Te_manual.Value;
        alpha_initial = 0.5;
        amp_initial = max(thomson_spectrum(fit_mask));

        num_starts = handles.edit_multistarts.Value;
        best_p_fit = []; best_resnorm = inf; best_jacobian = [];
        options = optimoptions('lsqcurvefit', 'Display', 'final'); % Display set to 'final'

        for start_idx = 1:num_starts
            if start_idx == 1
                p0_current = [Te_initial, alpha_initial, amp_initial];
            else
                p0_current = [Te_initial*(0.7 + 1.3*rand()), 2*rand(), amp_initial*(0.8 + 0.4*rand())];
            end

            lb = [0.1, 0.01, 0];
            ub = [15.0, 5.0, 3 * amp_initial];

            fit_fun = @(params, x) coherent_fitter_direct(params, x, laser_wl);

            try
                [p_fit_current, resnorm_current, ~, exitflag, ~, ~, jacobian_current] = lsqcurvefit(fit_fun, p0_current, wavelength_nm(fit_mask), thomson_spectrum(fit_mask), lb, ub, options);

                if exitflag > 0 && resnorm_current < best_resnorm
                    best_p_fit = p_fit_current;
                    best_resnorm = resnorm_current;
                    best_jacobian = jacobian_current;
                end
            catch ME
                fprintf('Multi-start %d failed: %s\n', start_idx, ME.message);
            end
        end

        if isempty(best_p_fit)
            uialert(ancestor(hObject, 'figure'), 'Shape fit failed to converge.', 'Error');
            update_status(handles, 'Shape fit failed', [1, 0, 0]);
            return;
        end

        Te_fit = best_p_fit(1);
        alpha_fit = best_p_fit(2);

        ne_fit = calculate_ne_from_alpha_Te(alpha_fit, Te_fit, laser_wl);

        % === NEW: Monte Carlo error estimation for Shape Fit ===
        fprintf('Calculating realistic parameter errors via Monte Carlo sampling...\n');

        % Define search ranges (±30% to capture typical variations)
        Te_range = linspace(Te_fit * 0.7, Te_fit * 1.3, 25);
        ne_range = linspace(ne_fit * 0.7, ne_fit * 1.3, 25);

        % Define acceptable SSR threshold (1.5x best fit)
        acceptable_ssr_threshold = best_resnorm * 1.5;
        fprintf('  SSR threshold for acceptable fits: %.2e (best fit SSR = %.2e)\n', acceptable_ssr_threshold, best_resnorm);

        acceptable_Te = [];
        acceptable_ne = [];
        amp_fit = best_p_fit(3);

        % Grid search over parameter space
        for Te_test = Te_range
            for ne_test = ne_range
                % Calculate alpha for this Te/ne combination
                alpha_test = calculate_alpha_from_ne(Te_test, ne_test, laser_wl);

                % Generate synthetic spectrum with these parameters
                test_spectrum = coherent_fitter_direct([Te_test, alpha_test, amp_fit], wavelength_nm, laser_wl);

                % Calculate unweighted residuals
                test_residuals = thomson_spectrum(fit_mask) - test_spectrum(fit_mask);
                test_ssr = sum(test_residuals.^2);

                % If this fit is within acceptable threshold, add to list
                if test_ssr < acceptable_ssr_threshold
                    acceptable_Te = [acceptable_Te, Te_test];
                    acceptable_ne = [acceptable_ne, ne_test];
                end
            end
        end

        % Calculate errors as half the range of acceptable values
        if ~isempty(acceptable_Te)
            Te_fit_err = (max(acceptable_Te) - min(acceptable_Te)) / 2;
            ne_fit_err = (max(acceptable_ne) - min(acceptable_ne)) / 2;
            fprintf('  Monte Carlo results: %d parameter sets within acceptable SSR\n', length(acceptable_Te));
            fprintf('  Te range: %.2f to %.2f eV (error = %.2f eV, %.1f%%)\n', min(acceptable_Te), max(acceptable_Te), Te_fit_err, 100*Te_fit_err/Te_fit);
            fprintf('  ne range: %.2e to %.2e m^-3 (error = %.2e m^-3, %.1f%%)\n', min(acceptable_ne), max(acceptable_ne), ne_fit_err, 100*ne_fit_err/ne_fit);
        else
            % Fallback: use 15% as empirically observed typical error
            Te_fit_err = Te_fit * 0.15;
            ne_fit_err = ne_fit * 0.15;
            fprintf('  Warning: No parameter sets found within threshold. Using fallback 15%% errors.\n');
        end

        fitted_spectrum = coherent_fitter_direct(best_p_fit, wavelength_nm, laser_wl);

        residuals = thomson_spectrum(fit_mask) - fitted_spectrum(fit_mask);
        ssr_shape = sum(residuals.^2);
        tss = sum((thomson_spectrum(fit_mask) - mean(thomson_spectrum(fit_mask))).^2);
        r_squared_shape = 1 - ssr_shape/tss;

        fprintf('Shape Fit Results: Te=%.2f +/- %.2f eV, alpha=%.2f, ne_inst=%.2e +/- %.2e m-3, SSR=%.2e, R²=%.3f\n', Te_fit, Te_fit_err, alpha_fit, ne_fit, ne_fit_err, ssr_shape, r_squared_shape);

        handles.shape_fit_ssr = ssr_shape;
        handles.shape_fit_r2 = r_squared_shape;
        handles.shape_fit_results.Te = Te_fit;
        handles.shape_fit_results.alpha = alpha_fit;
        handles.shape_fit_results.ne_inst = ne_fit;

        handles.edit_Te_manual.Value = Te_fit;
        handles.edit_alpha.Value = alpha_fit;
        handles.edit_ne_manual.Value = ne_fit;
        handles.edit_gamma.Value = 1.0;

        if isfield(handles, 'h_coherent_fit') && isvalid(handles.h_coherent_fit)
            set(handles.h_coherent_fit, 'YData', fitted_spectrum, 'DisplayName', 'Shape Fit', 'Visible', 'on');
            if isfield(handles, 'h_gauss_fit') && isvalid(handles.h_gauss_fit), set(handles.h_gauss_fit, 'Visible','off'); end
        else
            cla(handles.ax4);
            handles.h_data_points = plot(handles.ax4, wavelength_nm(fit_mask), thomson_spectrum(fit_mask), 'b-', 'DisplayName', 'Data');
            hold(handles.ax4, 'on');
            handles.h_coherent_fit = plot(handles.ax4, wavelength_nm, fitted_spectrum, 'g-', 'LineWidth', 2, 'DisplayName', 'Shape Fit');
            grid(handles.ax4, 'on'); legend(handles.ax4, 'Location', 'best');
            xlabel(handles.ax4, 'Wavelength (nm)'); ylabel(handles.ax4, 'Intensity (counts)');
        end

        title(handles.ax4, sprintf('Shape Fit (MC errors): Te=%.2f \\pm %.2f eV, n_{e,inst}=%.2e \\pm %.2e m^{-3}, \\alpha=%.2f\nSSR=%.2e, R^2=%.3f', Te_fit, Te_fit_err, ne_fit, ne_fit_err, alpha_fit, ssr_shape, r_squared_shape));

        update_status(handles, 'Ready', [0, 0.5, 0]);
        guidata(hObject, handles);
    end

    function scaled_spectrum = coherent_fitter_direct(p, wavelength_nm, laser_wl_nm)
        Te_eV = p(1);
        alpha = p(2);
        amplitude = p(3);

        theoretical_shape = calculate_coherent_ts(wavelength_nm, Te_eV, alpha, laser_wl_nm);
        max_theo = max(theoretical_shape);

        if max_theo > 1e-9
            normalized_shape = theoretical_shape / max_theo;
        else
            normalized_shape = zeros(size(wavelength_nm));
        end

        scaled_spectrum = amplitude * normalized_shape;
    end

    function shift_fit_left_callback(hObject, ~)
        handles = guidata(hObject);
        handles.fit_pixel_shift = handles.fit_pixel_shift - 1;
        apply_fit_shift(handles);
        guidata(hObject, handles);
    end

    function shift_fit_right_callback(hObject, ~)
        handles = guidata(hObject);
        handles.fit_pixel_shift = handles.fit_pixel_shift + 1;
        apply_fit_shift(handles);
        guidata(hObject, handles);
    end

    function apply_fit_shift(handles)
        if isfield(handles, 'h_coherent_fit') && isvalid(handles.h_coherent_fit) && strcmp(get(handles.h_coherent_fit, 'Visible'), 'on')
            fit_handle = handles.h_coherent_fit;
        elseif isfield(handles, 'h_gauss_fit') && isvalid(handles.h_gauss_fit)
            fit_handle = handles.h_gauss_fit;
        else
            return;
        end

        original_ydata = get(fit_handle, 'UserData');
        if isempty(original_ydata)
            original_ydata = get(fit_handle, 'YData');
            set(fit_handle, 'UserData', original_ydata);
        end

        shifted_ydata = circshift(original_ydata, handles.fit_pixel_shift);

        if handles.fit_pixel_shift > 0
            shifted_ydata(1:handles.fit_pixel_shift) = 0;
        elseif handles.fit_pixel_shift < 0
            shifted_ydata(end+handles.fit_pixel_shift+1:end) = 0;
        end

        set(fit_handle, 'YData', shifted_ydata);
        handles.label_shift.Text = sprintf('Shift: %d px', handles.fit_pixel_shift);
    end

    function reanalyze_all_callback(hObject,~)
        handles = guidata(hObject);
        update_status(handles, 'Performing intensity-calibrated Gaussian fit...', [1, 0.5, 0]);
        drawnow;
        disp('--- Intensity-Calibrated Gaussian Fit ---');

        % Get current parameters
        wavelength_nm = handles.plot_data.wavelength_nm;
        thomson_spectrum = handles.current_ts_spectrum;
        exclusion_cols = handles.edit_start.Value:handles.edit_end.Value;
        fit_mask = ~ismember(1:length(wavelength_nm), exclusion_cols);
        laser_wl = handles.params.laser_wavelength;

        % Fit a Gaussian to get FWHM and amplitude
        gauss_model = @(p, x) p(1) * exp(-4 * log(2) * ((x - p(2)).^2 / p(3)^2)) + p(4);

        % Multi-start Gaussian fitting
        num_starts = handles.edit_multistarts.Value;
        best_p_fit_gauss = [];
        best_resnorm_gauss = inf;
        best_jacobian = [];

        for start_idx = 1:num_starts
            if start_idx == 1
                p0_gauss = [max(thomson_spectrum(fit_mask)), laser_wl, 2.0, 0];
            else
                p0_gauss = [max(thomson_spectrum(fit_mask))*(0.8+0.4*rand()), laser_wl+0.2*randn(), 1.0+2.0*rand(), 0];
            end

            try
                [p_fit_current, resnorm_current, ~, exitflag, ~, ~, jacobian_current] = lsqcurvefit(gauss_model, p0_gauss, wavelength_nm(fit_mask), thomson_spectrum(fit_mask), [], [], optimset('Display','off'));

                if exitflag > 0 && resnorm_current < best_resnorm_gauss
                    best_p_fit_gauss = p_fit_current;
                    best_resnorm_gauss = resnorm_current;
                    best_jacobian = jacobian_current;
                end
            catch ME
                fprintf('  Multi-start %d failed: %s\n', start_idx, ME.message);
            end
        end

        if isempty(best_p_fit_gauss)
            uialert(ancestor(hObject, 'figure'), 'Gaussian fit failed to converge.', 'Error');
            update_status(handles, 'Fit failed', [1, 0, 0]);
            return;
        end

        p_fit_gauss = best_p_fit_gauss;

        % Extract fit parameters
        fwhm_fit = p_fit_gauss(3);
        amp_fit = p_fit_gauss(1);

        % Calculate Te from FWHM
        me_c2_eV = 0.511e6;
        Te_fit = me_c2_eV / (32 * log(2)) * (fwhm_fit / laser_wl)^2;

        % Calculate ne from intensity calibration
        gauss_fit_spectrum = gauss_model(p_fit_gauss, wavelength_nm);
        total_counts = sum(gauss_fit_spectrum);

        calib_factor = handles.plot_data.raman_results.calibration_factor;
        total_energy_thomson = (handles.edit_energy_t.Value/1000) * handles.edit_shots_t.Value;
        sigma_T_diff = 7.94e-30;

        ne_fit = total_counts / (calib_factor * total_energy_thomson * sigma_T_diff);

        % === NEW: Monte Carlo error estimation for Gaussian Fit ===
        fprintf('Calculating realistic parameter errors via Monte Carlo sampling...\n');

        % Define search ranges (±30% to capture typical variations)
        Te_range = linspace(Te_fit * 0.7, Te_fit * 1.3, 25);
        ne_range = linspace(ne_fit * 0.7, ne_fit * 1.3, 25);

        % Define acceptable SSR threshold (1.5x best fit)
        acceptable_ssr_threshold = best_resnorm_gauss * 1.5;
        fprintf('  SSR threshold for acceptable fits: %.2e (best fit SSR = %.2e)\n', acceptable_ssr_threshold, best_resnorm_gauss);

        acceptable_Te = [];
        acceptable_ne = [];
        center_wl = p_fit_gauss(2);
        baseline = p_fit_gauss(4);

        % Grid search over parameter space
        for Te_test = Te_range
            for ne_test = ne_range
                % Calculate FWHM from Te: Te = me_c2/(32*ln2) * (FWHM/laser_wl)^2
                % => FWHM = sqrt(Te * 32*ln2 / me_c2) * laser_wl
                fwhm_test = sqrt(Te_test * 32 * log(2) / me_c2_eV) * laser_wl;

                % Calculate amplitude from ne
                % ne = total_counts / (calib * energy * sigma)
                % total_counts = sum(A * exp(...)) ≈ A * sqrt(pi) * FWHM / sqrt(4*ln2)
                % So: A = ne * (calib * energy * sigma) / (sqrt(pi) * FWHM / sqrt(4*ln2))
                counts_needed = ne_test * calib_factor * total_energy_thomson * sigma_T_diff;
                gauss_integral_factor = sqrt(pi) * fwhm_test / sqrt(4 * log(2));
                amp_test = counts_needed / gauss_integral_factor;

                % Generate synthetic Gaussian spectrum
                test_spectrum = gauss_model([amp_test, center_wl, fwhm_test, baseline], wavelength_nm);

                % Calculate residuals
                test_residuals = thomson_spectrum(fit_mask) - test_spectrum(fit_mask);
                test_ssr = sum(test_residuals.^2);

                % If this fit is within acceptable threshold, add to list
                if test_ssr < acceptable_ssr_threshold
                    acceptable_Te = [acceptable_Te, Te_test];
                    acceptable_ne = [acceptable_ne, ne_test];
                end
            end
        end

        % Calculate errors as half the range of acceptable values
        if ~isempty(acceptable_Te)
            Te_fit_err = (max(acceptable_Te) - min(acceptable_Te)) / 2;
            ne_fit_err = (max(acceptable_ne) - min(acceptable_ne)) / 2;
            fprintf('  Monte Carlo results: %d parameter sets within acceptable SSR\n', length(acceptable_Te));
            fprintf('  Te range: %.2f to %.2f eV (error = %.2f eV, %.1f%%)\n', min(acceptable_Te), max(acceptable_Te), Te_fit_err, 100*Te_fit_err/Te_fit);
            fprintf('  ne range: %.2e to %.2e m^-3 (error = %.2e m^-3, %.1f%%)\n', min(acceptable_ne), max(acceptable_ne), ne_fit_err, 100*ne_fit_err/ne_fit);
        else
            % Fallback: use 15% as empirically observed typical error
            Te_fit_err = Te_fit * 0.15;
            ne_fit_err = ne_fit * 0.15;
            fprintf('  Warning: No parameter sets found within threshold. Using fallback 15%% errors.\n');
        end

        % Calculate R^2
        residuals = thomson_spectrum(fit_mask) - gauss_fit_spectrum(fit_mask);
        ssr = sum(residuals.^2);
        tss = sum((thomson_spectrum(fit_mask) - mean(thomson_spectrum(fit_mask))).^2);
        r_squared = 1 - ssr/tss;

        fprintf('Intensity-Calibrated Gaussian Fit Results:\n');
        fprintf('  Te = %.2f ± %.2f eV\n', Te_fit, Te_fit_err);
        fprintf('  ne = %.2e ± %.2e m^-3\n', ne_fit, ne_fit_err);
        fprintf('  FWHM = %.3f nm\n', fwhm_fit);
        fprintf('  SSR = %.2e, R² = %.3f\n', ssr, r_squared);

        % Update GUI fields
        handles.edit_Te_manual.Value = Te_fit;
        handles.edit_ne_manual.Value = ne_fit;
        alpha_new = calculate_alpha_from_ne(Te_fit, ne_fit, laser_wl);
        handles.edit_alpha.Value = alpha_new;
        handles.edit_gamma.Value = 1.0;

        % Update plot
        cla(handles.ax4);
        if strcmp(handles.plot_style, 'line')
            handles.h_data_points = plot(handles.ax4, wavelength_nm(fit_mask), thomson_spectrum(fit_mask), 'b-', 'DisplayName', 'Data');
        else
            handles.h_data_points = plot(handles.ax4, wavelength_nm(fit_mask), thomson_spectrum(fit_mask), 'bo', 'MarkerSize', 4, 'LineStyle', 'none', 'DisplayName', 'Data');
        end
        hold(handles.ax4, 'on');
        handles.h_gauss_fit = plot(handles.ax4, wavelength_nm, gauss_fit_spectrum, 'r-', 'LineWidth', 2, 'DisplayName', 'Intensity Cal Gaussian');
        grid(handles.ax4, 'on');
        legend(handles.ax4, 'Location', 'best');
        xlabel(handles.ax4, 'Wavelength (nm)');
        ylabel(handles.ax4, 'Intensity (counts)');
        title(handles.ax4, sprintf('Int. Cal. Gaussian (MC errors) (SSR=%.2e, R²=%.3f)\nTe = %.2f \\pm %.2f eV, ne = %.2e \\pm %.2e m^{-3}', ssr, r_squared, Te_fit, Te_fit_err, ne_fit, ne_fit_err));

        update_status(handles, 'Ready', [0, 0.5, 0]);
        guidata(hObject, handles);
        disp('--- Intensity-calibrated Gaussian fit complete. ---');
    end
    function toggle_plot_style_callback(hObject, ~)
        handles = guidata(hObject);
        if ~isfield(handles, 'h_data_points') || ~isvalid(handles.h_data_points), return; end
        if strcmp(handles.plot_style, 'marker')
            set(handles.h_data_points, 'Marker', 'none', 'LineStyle', '-');
            handles.plot_style = 'line';
        else
            set(handles.h_data_points, 'Marker', 'o', 'LineStyle', 'none');
            handles.plot_style = 'marker';
        end
        guidata(hObject, handles);
    end

    function recalculate_calib_callback(hObject, ~)
        handles = guidata(hObject);
        update_status(handles, 'Recalibrating from Raman (multi-start optimization)...', [1, 0.5, 0]);
        drawnow;

        pixels_exp = 1:size(handles.params.raman_img_raw, 2);

        % Stage 0: Intelligent peak-based calibration
        disp('=== Stage 0: Intelligent Peak-Based Wavelength Calibration ===');
        [a_improved, b_improved] = intelligent_peak_calibration_full(pixels_exp, handles.plot_data.raman_results.experimental_spectrum, {'N2', 'O2'}, handles.params.mixing_ratios, handles.params.pressure_Pa, handles.params.temp_K_initial, handles.params.linewidth_initial, handles.params.a_initial, handles.params.b_initial);

        handles.params.a_initial = a_improved;
        handles.params.b_initial = b_improved;
        handles.edit_a.Value = a_improved;
        handles.edit_b.Value = b_improved;

        % Stage 1: Multi-start optimization
        disp('=== Stage 1: Multi-start Physical Parameter Optimization ===');
        wavelength_exp_new = a_improved * pixels_exp + b_improved;
        total_energy_raman = handles.params.energy_J_raman * handles.params.shots_raman;
        theoretical_total_initial = calculate_theoretical_spectrum(wavelength_exp_new, {'N2', 'O2'}, handles.params.mixing_ratios, handles.params.pressure_Pa, handles.params.temp_K_initial, handles.params.linewidth_initial, handles.params.laser_wavelength);
        experimental_integral = trapz(wavelength_exp_new, handles.plot_data.raman_results.experimental_spectrum);
        theoretical_integral = trapz(wavelength_exp_new, theoretical_total_initial);
        if theoretical_integral > 0
            initial_calibration = experimental_integral / (total_energy_raman * theoretical_integral);
        else
            initial_calibration = 1e-15;
        end

        num_starts = handles.edit_multistarts.Value;
        best_params_raman = [];
        best_resnorm = inf;
        T_min = 273.15 + 17;
        T_max = 273.15 + 26;

        for start_idx = 1:num_starts
            if start_idx == 1
                initial_guess = [initial_calibration, handles.params.linewidth_initial, handles.params.temp_K_initial];
            else
                initial_guess = [initial_calibration * (0.6 + 0.8*rand()), handles.params.linewidth_initial * (0.7 + 0.8*rand()), T_min + (T_max - T_min) * rand()];
            end
            lower_bounds = [initial_calibration*0.4, handles.params.linewidth_initial*0.5, T_min];
            upper_bounds = [initial_calibration*1.6, handles.params.linewidth_initial*1.5, T_max];
            options = optimoptions('lsqnonlin', 'Display', 'final'); % Display set to 'final'
            try
                [current_fit_params, current_resnorm, ~, exitflag] = lsqnonlin(@(params) physical_objective_function_improved_full(params, wavelength_exp_new, handles.plot_data.raman_results.experimental_spectrum, {'N2', 'O2'}, handles.params.mixing_ratios, handles.params.pressure_Pa, total_energy_raman, handles.params.laser_wavelength), initial_guess, lower_bounds, upper_bounds, options);
                if current_resnorm < best_resnorm && exitflag > 0
                    best_params_raman = current_fit_params;
                    best_resnorm = current_resnorm;
                end
            catch ME
                fprintf('Start %d failed: %s\n', start_idx, ME.message);
            end
        end

        if ~isempty(best_params_raman)
            handles.params.linewidth_initial = best_params_raman(2);
            handles.params.temp_K_initial = best_params_raman(3);
            handles.edit_temp_K.Value = best_params_raman(3);
            handles.edit_linewidth.Value = best_params_raman(2);
        end

        current_params = handles.params;
        current_params.temp_K_initial = handles.edit_temp_K.Value;
        current_params.mixing_ratios = [handles.edit_n2_frac.Value, handles.edit_o2_frac.Value];

        % UPDATE params with current GUI values BEFORE calling run_full_analysis
        handles.params.shots_raman = handles.edit_shots_r.Value;
        handles.params.energy_J_raman = handles.edit_energy_r.Value / 1000;
        handles.params.shots_thomson = handles.edit_shots_t.Value;
        handles.params.energy_J_thomson = handles.edit_energy_t.Value / 1000;

        [results, plot_data] = run_full_analysis(handles.params);

        handles.plot_data = plot_data;
        handles.current_ts_spectrum = plot_data.thomson_spectrum_processed;
        handles.current_results.uncorrected = results.uncorrected;
        handles.current_deconv = results.deconv_results;

        % UPDATE GUI FIELDS
        handles.edit_Te_manual.Value = plot_data.deconv_results.Te_eV_corrected;
        handles.edit_ne_manual.Value = plot_data.deconv_results.ne_corrected;
        alpha_new = calculate_alpha_from_ne(plot_data.deconv_results.Te_eV_corrected, plot_data.deconv_results.ne_corrected, handles.params.laser_wavelength);
        handles.edit_alpha.Value = alpha_new;

        update_plots(handles);
        update_status(handles, 'Ready', [0, 0.5, 0]);
        guidata(hObject, handles);
    end

    function update_plots(handles)
        set(handles.h_thomson_img, 'CData', medfilt2(handles.plot_data.thomson_img_processed, [handles.plot_data.mfp handles.plot_data.mfp]));
        set(handles.h_raman_img, 'CData', handles.plot_data.raman_img_processed);

        set(handles.h_raman_spec, 'XData', handles.plot_data.raman_results.wavelength_axis, 'YData', handles.plot_data.raman_results.experimental_spectrum);
        model_scaled = handles.plot_data.raman_results.theoretical_spectrum_model * (max(handles.plot_data.raman_results.experimental_spectrum)/max(handles.plot_data.raman_results.theoretical_spectrum_model));
        set(handles.h_raman_model, 'XData', handles.plot_data.raman_results.wavelength_axis, 'YData', model_scaled);
        set(handles.h_calib_title, 'String', sprintf('Raman Spectrum (Calib. Factor: %.2e)', handles.plot_data.raman_results.calibration_factor));

        fit_mask = ~ismember(1:length(handles.plot_data.wavelength_nm), handles.params.exclusion_cols);
        gauss_model = @(p, x) p(1) * exp(-4 * log(2) * ((x - p(2)).^2 / p(3)^2)) + p(4);
        set(handles.h_data_points, 'XData', handles.plot_data.wavelength_nm(fit_mask), 'YData', handles.plot_data.thomson_spectrum_processed(fit_mask));
        if isfield(handles, 'h_gauss_fit') && isvalid(handles.h_gauss_fit)
            set(handles.h_gauss_fit, 'XData', handles.plot_data.wavelength_nm, 'YData', gauss_model(handles.plot_data.gauss_fit_params, handles.plot_data.wavelength_nm));
        end
        set(handles.h_title, 'String', sprintf('Direct TS Gaussian Fit (R^2 = %.3f)\nTe = %.2f \\pm %.2f eV, ne = %.2e \\pm %.2e m^{-3}', handles.plot_data.r_squared_gauss, handles.current_deconv.Te_eV_corrected, handles.current_deconv.Te_eV_error, handles.current_deconv.ne_corrected, handles.current_deconv.ne_error));
    end

    function decouple_plots_callback(hObject,~)
        handles = guidata(hObject);
        disp('Decoupling plots into separate figures...');
        figure('Name', 'Thomson Image'); new_ax = axes(); copyobj(allchild(handles.ax1), new_ax); title(new_ax,'Thomson Image'); axis(new_ax, 'image'); colormap(new_ax, 'jet');
        figure('Name', 'Raman Image'); new_ax = axes(); copyobj(allchild(handles.ax2), new_ax); title(new_ax,'Raman Image'); axis(new_ax, 'image'); colormap(new_ax, 'jet');
        figure('Name', 'Raman Spectrum'); new_ax = axes(); copyobj(allchild(handles.ax3), new_ax); title(new_ax,get(get(handles.ax3, 'Title'), 'String')); grid(new_ax, 'on'); legend(new_ax);
        figure('Name', 'Thomson Spectrum'); new_ax = axes(); copyobj(allchild(handles.ax4), new_ax); title(new_ax,get(get(handles.ax4, 'Title'), 'String')); grid(new_ax, 'on'); legend(new_ax);
    end

    function plot_eedf_callback(hObject,~)
        handles = guidata(hObject);
        disp('Plotting EEDF with current data...');
        generate_eedf_plot(handles.plot_data.wavelength_nm, handles.current_ts_spectrum, handles.current_results.uncorrected.ne_m3, 532);
    end

    function save_ts_data_callback(hObject,~)
        handles = guidata(hObject);
        update_status(handles, 'Saving TS data...', [0, 0, 1]);

        try
            % 1. Determine which fit is currently displayed
            fit_handle = [];
            fit_type_raw = '';
            ssr = NaN;

            % Check if the coherent fit (Shape, Intensity, Recalc, etc.) is visible
            if isfield(handles, 'h_coherent_fit') && isvalid(handles.h_coherent_fit) && strcmp(get(handles.h_coherent_fit, 'Visible'), 'on')
                fit_handle = handles.h_coherent_fit;
                fit_type_raw = get(fit_handle, 'DisplayName');

                % Determine which type of coherent fit based on DisplayName
                if contains(fit_type_raw, 'Intensity Cal')
                    if isfield(handles, 'intensity_fit_ssr')
                        ssr = handles.intensity_fit_ssr;
                    end
                    fit_type_raw = 'Intensity Cal';
                elseif contains(fit_type_raw, 'Shape Fit')
                    if isfield(handles, 'shape_fit_ssr')
                        ssr = handles.shape_fit_ssr;
                    end
                    fit_type_raw = 'Shape Fit';
                elseif contains(fit_type_raw, 'Recalculated')
                    fit_type_raw = 'Recalculated';
                    ssr = NaN; % Recalc doesn't have SSR
                else
                    fit_type_raw = 'Coherent Fit';
                end

                % If not, check if the Gaussian fit is visible
            elseif isfield(handles, 'h_gauss_fit') && isvalid(handles.h_gauss_fit) && strcmp(get(handles.h_gauss_fit, 'Visible'), 'on')
                fit_handle = handles.h_gauss_fit;
                fit_type_raw = 'Gaussian';
                ssr = handles.plot_data.ssr_gauss;
            end

            if isempty(fit_handle)
                uialert(ancestor(hObject,'figure'), 'No active fit found to save.', 'Save Error');
                update_status(handles, 'Save failed: No active fit.', [1, 0, 0]);
                return;
            end

            % 2. Get CURRENT data for the file content and name
            fitted_spectrum = get(fit_handle, 'YData');
            Te = handles.edit_Te_manual.Value;  % Use current GUI values
            ne = handles.edit_ne_manual.Value;  % Use current GUI values
            alpha = handles.edit_alpha.Value;   % Get alpha too

            % Clean up the fit type string for the filename
            fit_type_clean = strrep(fit_type_raw, ' ', ''); % Removes spaces

            % 3. Construct the filename and the full save path
            [save_path, ~, ~] = fileparts(handles.thomson_file_original);

            if isfinite(ssr)
                file_name = sprintf('TS_%s_Te%.2feV_ne%.2e_alpha%.2f_SSR%.2e.csv', ...
                    fit_type_clean, Te, ne, alpha, ssr);
            else
                file_name = sprintf('TS_%s_Te%.2feV_ne%.2e_alpha%.2f.csv', ...
                    fit_type_clean, Te, ne, alpha);
            end

            full_file_path = fullfile(save_path, file_name);

            % 4. Prepare the data in a table for saving
            wavelength_nm = handles.plot_data.wavelength_nm(:);
            experimental_spectrum = handles.current_ts_spectrum(:);
            fitted_spectrum = fitted_spectrum(:);

            data_table = table(wavelength_nm, experimental_spectrum, fitted_spectrum);
            data_table.Properties.VariableNames = {'Wavelength_nm', 'Experimental_Intensity', 'Fitted_Intensity'};

            % 5. Write the table to a CSV file
            writetable(data_table, full_file_path);

            % 6. Notify the user of success
            fprintf('Successfully saved TS data to: %s\n', full_file_path);
            uialert(ancestor(hObject,'figure'), sprintf('Data saved successfully to:\n%s', full_file_path), 'Save Successful');
            update_status(handles, 'Ready', [0, 0.5, 0]);

        catch ME
            fprintf('Error saving TS data: %s\n', ME.message);
            uialert(ancestor(hObject,'figure'), ['An error occurred while saving: ' ME.message], 'Save Error');
            update_status(handles, 'Save failed.', [1, 0, 0]);
        end
    end
end
% =========================================================================
%               HELPER, ANALYSIS & PHYSICS FUNCTIONS
% =========================================================================

function fitted_spectrum = intensity_calibration_objective(params, wavelength_nm, laser_wl_nm, calib_factor, total_energy, sigma_diff)
Te_eV = params(1);
ne = params(2);

alpha = calculate_alpha_from_ne(Te_eV, ne, laser_wl_nm);
shape = calculate_coherent_ts(wavelength_nm, Te_eV, alpha, laser_wl_nm);

sum_shape = sum(shape);
if sum_shape > 1e-12
    normalized_shape = shape / sum_shape;
else
    normalized_shape = zeros(size(wavelength_nm));
end

total_signal = calib_factor * total_energy * sigma_diff * ne;
fitted_spectrum = total_signal * normalized_shape;
end

function ts_spectrum = calculate_coherent_ts(wavelength_nm, Te_eV, alpha, laser_wl_nm)
kB_J = 1.38065e-23; kB_eV = 8.617343e-5; c = 299792458; m_e = 9.1094e-31;
angle_rad = 90 * (pi/180); laser_wl_m = laser_wl_nm * 1e-9;
k_wave = (4*pi / laser_wl_m) * sin(angle_rad / 2);
ve = sqrt(2 * kB_J * (Te_eV / kB_eV) / m_e);
w_laser = 2*pi*c / laser_wl_m; w_axis = 2*pi*c ./ (wavelength_nm * 1e-9);
w_shift = w_axis - w_laser; xe = w_shift ./ (k_wave * ve);
W = 1 - 2.*xe.*dawson(xe) - 1i.*sqrt(pi).*xe.*exp(-xe.^2);
denom = (1 + alpha.^2.*W) .* conj(1 + alpha.^2.*W);
Se = exp(-xe.^2) ./ (sqrt(pi) * denom);
ts_spectrum = real(Se); ts_spectrum(isnan(ts_spectrum)) = 0;
end

function alpha = calculate_alpha_from_ne(Te_eV, ne, laser_wl_nm)
q_e = 1.6022e-19; eps0 = 8.8541878128e-12;
k_wave = (4 * pi / (laser_wl_nm * 1e-9)) * sin((90 * pi / 180) / 2);
kT_J = Te_eV * q_e;
lambda_D = sqrt(eps0 * kT_J / (ne * q_e^2));
alpha = 1 / (k_wave * lambda_D);
end

function ne = calculate_ne_from_alpha_Te(alpha, Te_eV, laser_wl_nm)
q_e = 1.6022e-19; eps0 = 8.8541878128e-12; kB_J = 1.38065e-23; kB_eV = 8.617343e-5;
k_wave = (4*pi / (laser_wl_nm*1e-9)) * sin((90*pi/180) / 2);
ne = alpha^2 * k_wave^2 * (eps0 * kB_J * (Te_eV/kB_eV)) / q_e^2;
end

function deconv_results = amplitude_fit_deconvolution(wl, data, calib, lw, params, fwhm0, amp0)
laser_wl = 532; meV = 0.511e6; sigma_T = 7.94e-30;
opts = optimset('Display', 'final', 'TolFun', 1e-9, 'TolX', 1e-9); % Display set to 'final'
wing_mask = abs(wl - laser_wl) > 1.5;
wing_model = @(p, x) p(1) * exp(-4 * log(2) * ((x - laser_wl) ./ p(2)).^2);

% Get Jacobian from the wing fit to calculate error on FWHM, which gives Te error
[p_wing, resnorm_wing, ~, ~, ~, ~, jacobian_wing] = lsqcurvefit(wing_model, [max(data(wing_mask)), fwhm0], wl(wing_mask), data(wing_mask), [0, 0.1], [amp0 * 2, 10], opts);
fwhm_fixed = p_wing(2);
Te_wings = meV / (32 * log(2)) * (fwhm_fixed / laser_wl)^2;

% Calculate Te_err from fwhm_err
Te_err = 0; % Default value
dof_wing = sum(wing_mask) - length(p_wing);
if dof_wing > 0 && rcond(full(jacobian_wing' * jacobian_wing)) > 1e-15
    covB_wing = resnorm_wing / dof_wing * inv(jacobian_wing' * jacobian_wing);
    p_wing_err = full(sqrt(diag(covB_wing)));
    fwhm_err = p_wing_err(2);
    % Propagate error: Te = C*fwhm^2  => dTe/Te = 2*dfwhm/fwhm
    Te_err = Te_wings * (2 * fwhm_err / fwhm_fixed);
end

p0 = [amp0, amp0 * 0.1, 500]; lb = [0, 0, 300]; ub = [amp0 * 3, amp0 * 2, 1500];
obj_func = @(p) objective_func_deconv(p, wl, data, fwhm_fixed, lw, laser_wl);
[p_fit, resnorm, ~, ~, ~, ~, jacobian] = lsqnonlin(obj_func, p0, lb, ub, opts);

% Calculate ne_err from amplitude error
ne_err = 0; % Default value
dof = length(data) - length(p_fit);
if dof > 0 && rcond(full(jacobian' * jacobian)) > 1e-15
    covB = resnorm / dof * inv(jacobian' * jacobian);
    p_err = full(sqrt(diag(covB)));
    ne_err_fraction = p_err(1) / p_fit(1); % Error from the TS amplitude fit parameter
else
    ne_err_fraction = inf;
end

ts_amp = p_fit(1); rs_amp = p_fit(2); tg_fit = p_fit(3);
total_energy_thomson = params.energy_J_thomson * params.shots_thomson;

ts_component_fit = ts_amp * exp(-4*log(2)*((wl - laser_wl)./fwhm_fixed).^2);
total_counts_ts_component = sum(ts_component_fit);
ne = total_counts_ts_component / (calib * total_energy_thomson * sigma_T);
ne_err = ne * ne_err_fraction; % Apply fractional error to get absolute error

rs100 = calculateraman(wl, params.laser_wavelength, 'N2', tg_fit, params.pressure_Pa * 0.79, 'gaussian', lw) + ...
    calculateraman(wl, params.laser_wavelength, 'O2', tg_fit, params.pressure_Pa * 0.21, 'gaussian', lw);
peak100 = calib * total_energy_thomson * max(rs100);
air_frac = rs_amp / peak100;
[~, ts_comp, rs_comp] = objective_func_deconv(p_fit, wl, data, fwhm_fixed, lw, laser_wl);
combo = ts_comp + rs_comp;
r2 = 1 - sum((data - combo).^2) / sum((data - mean(data)).^2);
deconv_results = struct('Te_eV_corrected', Te_wings, 'ne_corrected', ne, ...
    'air_fraction_percent', air_frac * 100, 'gas_temperature_K', tg_fit, ...
    'thomson_component', ts_comp, 'raman_component', rs_comp, ...
    'combined_spectrum', combo, 'r_squared_deconv', r2, 'params', params, ...
    'Te_eV_error', Te_err, 'ne_error', ne_err);
end

function generate_eedf_plot(wavelength_nm, thomson_spectrum, ne, laser_wavelength)
figure('Name', 'Electron Energy Distribution Function (EEDF)');
c=299792458; m_e=9.1094e-31; q_e=1.6022e-19;
delta_lambda = wavelength_nm - laser_wavelength;
vk = (delta_lambda*1e-9*c) / (2*(laser_wavelength*1e-9)*sin(pi/2));
energy_eV = 0.5*m_e*vk.^2/q_e;
normalized_spectrum = thomson_spectrum / sum(thomson_spectrum);
blue_mask = wavelength_nm < laser_wavelength & thomson_spectrum > 0;
red_mask = wavelength_nm > laser_wavelength & thomson_spectrum > 0;
plot(energy_eV(blue_mask), normalized_spectrum(blue_mask), 'bo', 'DisplayName','Blue Wing'); hold on;
plot(energy_eV(red_mask), normalized_spectrum(red_mask), 'ro', 'DisplayName','Red Wing');
title(sprintf('EEDF (from uncorrected data, ne=%.2e m^{-3})', ne));
xlabel('Electron Energy (eV)'); ylabel('EEDF (normalized)');
grid on; set(gca, 'YScale', 'log'); legend('show'); ylim([1e-5, 1]);
end

function image_out = block_central_pixels(image_in, rows, cols)
roi=image_in(rows,:); n=size(roi,2); w=10;
baseline=linspace(mean(mean(roi(:,1:w))), mean(mean(roi(:,end-w+1:end))), n);
image_out=image_in;
image_out(:,cols)=repmat(baseline(cols), size(image_in,1), 1);
end

function spectrum_out = subtract_wing_baseline(spectrum_in)
n=length(spectrum_in); w=100;
baseline=linspace(mean(spectrum_in(1:w)), mean(spectrum_in(end-w+1:end)), n);
spectrum_out=max(spectrum_in-baseline, 0);
end

function theoretical_spectrum = calculate_theoretical_spectrum(wl, sp, r, P, T, lw,laser_nm)
theoretical_spectrum=zeros(size(wl));
for i=1:length(sp)
    theoretical_spectrum=theoretical_spectrum + calculateraman(wl,laser_nm, sp{i}, T, P*r(i), 'gaussian', lw);
end
end

function [resid, ts_comp, rs_comp] = objective_func_deconv(p, wl, data, fwhm, lw, laser_wl)
ts_amp=p(1); rs_amp=p(2); tg=p(3);
ts_comp=ts_amp*exp(-4*log(2)*((wl-laser_wl)./fwhm).^2);
rs_shape=calculateraman(wl, laser_wl, 'N2',tg,101325*0.79,'gaussian',lw)+calculateraman(wl, laser_wl, 'O2',tg,101325*0.21,'gaussian',lw);
if max(rs_shape)>0, rs_comp=rs_amp*rs_shape/max(rs_shape); else, rs_comp=zeros(size(wl)); end
resid=data-(ts_comp+rs_comp);
end

function residual = physical_objective_function_improved_full(params, wavelength, experimental_spectrum, species_list, mixing_ratios, P_atm, total_energy, laser_nm)
calibration_factor = params(1); linewidth = params(2); temperature = params(3);
if any(params <= 0), residual = ones(size(experimental_spectrum)) * 1e10; return; end
try
    theoretical_spectrum = calculate_theoretical_spectrum(wavelength, species_list, mixing_ratios, P_atm, temperature, linewidth,laser_nm);
    fitted_spectrum = calibration_factor * total_energy * theoretical_spectrum;
    residual = experimental_spectrum - fitted_spectrum;
catch, residual = ones(size(experimental_spectrum)) * 1e10; end
residual(~isfinite(residual)) = 1e10;
end

function [improved_a, improved_b] = intelligent_peak_calibration_full(pixels, experimental_spectrum, species_list, mixing_ratios, P_atm, T_ref, linewidth_ref, a_initial, b_initial)
fprintf('Starting intelligent peak-based calibration...\n');
exp_smooth = smoothdata(experimental_spectrum, 'gaussian', 5);
thresholds = [0.2, 0.3, 0.4, 0.5]; all_exp_peaks = [];
for thresh = thresholds
    [peaks, locs] = findpeaks(exp_smooth, 'MinPeakHeight', max(exp_smooth)*thresh, 'MinPeakDistance', 8, 'MinPeakProminence', max(exp_smooth)*thresh*0.3);
    if length(locs) >= 3, all_exp_peaks = [all_exp_peaks; [locs', peaks']]; end
end
if isempty(all_exp_peaks), fprintf('Warning: No significant peaks found. Using original calibration.\n'); improved_a = a_initial; improved_b = b_initial; return; end
unique_peaks = [];
for i = 1:size(all_exp_peaks, 1)
    if isempty(unique_peaks) || all(abs(all_exp_peaks(i,1) - unique_peaks(:,1)) > 5), unique_peaks = [unique_peaks; all_exp_peaks(i,:)]; end
end
exp_peak_pixels = unique_peaks(:,1); exp_peak_intensities = unique_peaks(:,2);
fprintf('Found %d experimental peaks\n', length(exp_peak_pixels));

wavelength_ref = a_initial * pixels + b_initial;
theoretical_ref = calculate_theoretical_spectrum(wavelength_ref, species_list, mixing_ratios, P_atm, T_ref, linewidth_ref, 532.0);
theo_smooth = smoothdata(theoretical_ref, 'gaussian', 5);
[theo_peaks, theo_locs] = findpeaks(theo_smooth, 'MinPeakHeight', max(theo_smooth)*0.2, 'MinPeakDistance', 8, 'MinPeakProminence', max(theo_smooth)*0.1);
if length(theo_locs) < 3, fprintf('Warning: Insufficient theoretical peaks found. Using original calibration.\n'); improved_a = a_initial; improved_b = b_initial; return; end
theo_peak_wavelengths = wavelength_ref(theo_locs);
fprintf('Found %d theoretical peaks\n', length(theo_peak_wavelengths));

best_a = a_initial; best_b = b_initial; best_score = inf;
[~, exp_strongest_idx] = sort(exp_peak_intensities, 'descend');
[~, theo_strongest_idx] = sort(theo_peaks, 'descend');
n_match = min([5, length(exp_strongest_idx), length(theo_locs)]);

if n_match >= 2
    for i = 1:min(3, n_match-1)
        for j = i+1:min(i+2, n_match)
            exp_pix1 = exp_peak_pixels(exp_strongest_idx(i)); exp_pix2 = exp_peak_pixels(exp_strongest_idx(j));
            theo_wav1 = theo_peak_wavelengths(theo_strongest_idx(i)); theo_wav2 = theo_peak_wavelengths(theo_strongest_idx(j));
            if abs(exp_pix2 - exp_pix1) > 10
                a_test = (theo_wav2 - theo_wav1) / (exp_pix2 - exp_pix1);
                b_test = theo_wav1 - a_test * exp_pix1;
                if a_test > 0.005 && a_test < 0.02 && abs(b_test) < 600
                    score = evaluate_calibration_quality_full(pixels, experimental_spectrum, species_list, mixing_ratios, P_atm, T_ref, linewidth_ref, a_test, b_test);
                    if score < best_score, best_a = a_test; best_b = b_test; best_score = score;
                        fprintf('  New best from strong peaks %d,%d: a=%.6f, b=%.4f, score=%.4f\n', i, j, a_test, b_test, score);
                    end
                end
            end
        end
    end
end

if best_score < inf
    a_range = linspace(best_a*0.9, best_a*1.1, 15); b_range = linspace(best_b-2, best_b+2, 15);
    for a_test = a_range
        for b_test = b_range
            score = evaluate_calibration_quality_full(pixels, experimental_spectrum, species_list, mixing_ratios, P_atm, T_ref, linewidth_ref, a_test, b_test);
            if score < best_score, best_a = a_test; best_b = b_test; best_score = score; end
        end
    end
end

wavelength_test = best_a * pixels + best_b;
theoretical_test = calculate_theoretical_spectrum(wavelength_test, species_list, mixing_ratios, P_atm, T_ref, linewidth_ref, 532.0);
exp_norm = experimental_spectrum / max(experimental_spectrum); theo_norm = theoretical_test / max(theoretical_test);
shifts = -5:0.5:5; best_corr = -inf; best_shift = 0;
for shift = shifts
    shifted_wavelength = wavelength_test + shift;
    shifted_theo = interp1(wavelength_test, theo_norm, shifted_wavelength, 'linear', 0);
    if max(shifted_theo) > 0
        corr_matrix = corrcoef(exp_norm, shifted_theo);
        if size(corr_matrix, 1) > 1 && corr_matrix(1,2) > best_corr, best_corr = corr_matrix(1,2); best_shift = shift; end
    end
end
if abs(best_shift) > 0.1, best_b = best_b + best_shift; fprintf('  Applied wavelength shift: %.3f nm, correlation improved to %.4f\n', best_shift, best_corr); end
improved_a = best_a; improved_b = best_b;
fprintf('Intelligent calibration completed: a=%.6f, b=%.4f\n', improved_a, improved_b);
end

function score = evaluate_calibration_quality_full(pixels, experimental_spectrum, species_list, mixing_ratios, P_atm, T_ref, linewidth_ref, a_test, b_test)
try
    wavelength_test = a_test * pixels + b_test;
    theoretical_test = calculate_theoretical_spectrum(wavelength_test, species_list, mixing_ratios, P_atm, T_ref, linewidth_ref, 532.0);
    if max(theoretical_test) <= 0, score = 1000; return; end
    exp_norm = experimental_spectrum / max(experimental_spectrum); theo_norm = theoretical_test / max(theoretical_test);
    corr_matrix = corrcoef(exp_norm, theo_norm);
    if size(corr_matrix, 1) > 1, correlation = corr_matrix(1,2); else, correlation = 0; end
    rms_diff = sqrt(mean((exp_norm - theo_norm).^2));
    exp_smooth = smoothdata(exp_norm, 'gaussian', 3); theo_smooth = smoothdata(theo_norm, 'gaussian', 3);
    [~, exp_locs] = findpeaks(exp_smooth, 'MinPeakHeight', 0.25, 'MinPeakDistance', 8);
    [~, theo_locs] = findpeaks(theo_smooth, 'MinPeakHeight', 0.25, 'MinPeakDistance', 8);
    peak_penalty = 0;
    if length(exp_locs) >= 2 && length(theo_locs) >= 2
        exp_wavelengths = wavelength_test(exp_locs); theo_wavelengths = wavelength_test(theo_locs);
        total_distance = 0; matched_peaks = 0;
        for i = 1:length(exp_wavelengths)
            [min_dist, ~] = min(abs(exp_wavelengths(i) - theo_wavelengths));
            if min_dist < 1.5, total_distance = total_distance + min_dist; matched_peaks = matched_peaks + 1; end
        end
        if matched_peaks > 0, peak_penalty = total_distance / matched_peaks; else, peak_penalty = 10; end
    else, peak_penalty = 5;
    end
    score = 0.5 * (1 - correlation) + 0.3 * rms_diff + 0.2 * peak_penalty;
catch
    score = 1000;
end
end

function ts_filtered = apply_ts_filters(ts_spectrum, remove_outliers, apply_medfilt, medfilt_window)
% Apply outlier removal and/or median filtering to TS spectrum
ts_filtered = ts_spectrum;

if remove_outliers
    ts_filtered = remove_outliers_iqr(ts_filtered);
end

if apply_medfilt
    % Ensure window is odd
    if mod(medfilt_window, 2) == 0
        medfilt_window = medfilt_window + 1;
    end
    ts_filtered = medfilt1(ts_filtered, medfilt_window);
end
end

function data_clean = remove_outliers_iqr(data)
% Remove extreme outliers using IQR method (3*IQR threshold)
q1 = quantile(data, 0.25);
q3 = quantile(data, 0.75);
iqr = q3 - q1;
lower_bound = q1 - 3 * iqr;
upper_bound = q3 + 3 * iqr;

% Replace outliers with NaN
data_clean = data;
outlier_mask = (data < lower_bound) | (data > upper_bound);
data_clean(outlier_mask) = NaN;

% Interpolate NaN values
if any(outlier_mask)
    valid_idx = find(~isnan(data_clean));
    if length(valid_idx) > 1
        data_clean = interp1(valid_idx, data_clean(valid_idx), 1:length(data_clean), 'linear', 'extrap');
    end
end
end

   