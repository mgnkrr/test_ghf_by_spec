%% ================================================================
%  GHF model evaluation with uncertainty (RAW vs THK-ADJ)
%  — Connected sinks, analytic expected SPEC/REC + sink-heneCIs —
%
%  Differences vs. prior:
%    • Orientation sanity: flip rasters if Y grid is descending
%    • Uncertainty: Expected SPEC/REC using Pr(GΔ>0) with σ from model
%      grids (or fallbacks) PLUS per-pixel σ(Gmin) from finite-difference.
%    • CIs from a sink-component bootstrap (resample connected sinks).
%    • Shaded CI bands (50%/90%).
%    • Sinks restrict ROI **with connectivity grouping**.
%
%  Inputs expected (from preprocessing):
%    datasets_for_gmin/GHF_*_interp.mat              (mean fields)
%    datasets_for_gmin/UNC_*_interp.mat              (σ fields, where applicable)
%    datasets_for_gmin/BMIN_Losing_interp.mat
%    datasets_for_gmin/BMAX_Losing_interp.mat
%    datasets_for_gmin/specularity.mat               (Q)
%    datasets_for_gmin/coldex_icethk.mat             (Xgrid, Ygrid, H)
%    datasets_for_gmin/Ts_interp.mat                 (Ts_interp)
%    datasets_for_gmin/Mb_interp.mat                 (Mb_interp)
%    datasets_for_gmin/mouginot_icevel.mat           (icevel)
%    datasets_for_gmin/sink_mask_new.mat (logical)   var: sink_mask
%  ================================================================
function S = bootstrap_ghf_component(cfg_overrides)

clearvars;
check_dependencies(struct('autoInstall',true,'needAMT',false,'verbose',true));

%% -------------------- CONFIG --------------------
cfg = struct( ...
  'spec_thresh',        0.20, ...
  'v_keep',             3, ...
  'target_pos_pct',     10, ...
  'outdir',             'figs_out', ...
  'overwrite',          true, ...
  'overlay_alpha',      0.75, ...
  'H_contour_interval', 250, ...
  'y_right',            true, ...
  'profile',            false, ...
  'to_single',          true, ...
  'keep_XY',            false, ...
  'debug',              true, ...
  'close_figs',         false  ...
);

% Flat GHF sweep
cfg.flat = struct( ...
  'enable',         true, ...
  'values',         25:1:90, ...
  'interp_enable',  true, ...
  'interp_points',  200, ...
  'interp_method',  'spline' ...
);

% Sinks (used to define ROI; WITH connectivity grouping)
cfg.sinks = struct( ...
  'enable',        true, ...
  'path',          'datasets_for_gmin/sink_mask_new.mat', ...
  'varname',       'sink_mask', ...
  'marker',        '.', ...
  'size',          4, ...
  'color',         [0.05 0.05 0.05], ...
  'alpha',         0.85, ...
  'random_thin',   1.0, ...
  'connectivity',  8 ...
);

% Uncertainty — analytic expectation + sink bootstrap
cfg.uncertainty = struct( ...
  'mode',        'analytic', ...  % (reserved)
  'bootstrap_mode', 'component', ... % 'pixel' | 'component'
  'n_boot',      500, ...         % # component resamples
  'band_main',   0.95, ...        % outer shaded band (e.g., 90%)
  'band_inner',  0.50 ...         % inner shaded band (e.g., 50%)
);

% Gmin uncertainty via finite-difference propagation (per pixel)
cfg.gmin_unc = struct( ...
  'enable',     true, ... 
  'sigma_Ts',   1.5, ...   % K van Den Broeke 2004
  'sigma_Mb',   0.02, ...  % Mb van de Berg
  'sigma_H',    50.0, ...  % m Duncan
  'delta_Ts',   0.5, ...
  'delta_Mb',   0.01, ...
  'delta_H',    5.0 ...
);

% Stage toggles
run = struct('load_data',true,'build_labels',true,'process_ghf',true,'evaluate',true);

% Plot toggles
cfg.plots = struct( ...
  'prrec_expected_adj',    true, ...
  'prrec_expected_raw',    false, ...
  'maps',                  true, ...
  'flat_markers',          false, ...
  'flat_fill',             true, ...
  'flat_fill_alpha',       0.10 ...
);

% Two-stat colored-CI figures to render (each cell is {Xstat, Ystat}) --> SPEC','REC','PREC','ACC','F1'
cfg.two_stat_pairs = { {'SPEC','REC'}};

% Global fallback per-model σ (mW/m^2) used only when no grid is present
fallback_err = struct('FoxMaule',21,'Hazzard',8,'Martos',10,'Shen',9,'Stal',12,'An',45,'Losing',25);

if nargin>=1 && ~isempty(cfg_overrides), cfg = override_struct(cfg, cfg_overrides); end
if ~exist(cfg.outdir,'dir'), mkdir(cfg.outdir); end
if cfg.profile, profile on; end

% Replace BOTH debug-normalization blocks with this single one:
if ~isfield(cfg,'debug') || (~isstruct(cfg.debug) && ~islogical(cfg.debug))
    cfg.debug = false;
end
if islogical(cfg.debug)
    cfg.debug = struct('enable',logical(cfg.debug), 'level',2, 'sample_k',10, ...
        'save_snapshots',false, 'snapshot_dir',fullfile(cfg.outdir,'debug'), ...
        'subset_roi',1e6, 'quicklook_plots',true, 'pause_on_fig',false, ...
        'progress_every',25,'quicklook_all',true,'quicklook_pause',false);
else
    def = struct('enable',false,'level',1,'sample_k',10,'save_snapshots',false, ...
        'snapshot_dir',fullfile(cfg.outdir,'debug'),'subset_roi',1e6, ...
        'quicklook_plots',false,'pause_on_fig',false,'progress_every',25, ...
        'quicklook_all',true,'quicklook_pause',false);
    fn = fieldnames(def); for ii=1:numel(fn), if ~isfield(cfg.debug,fn{ii}), ...
                cfg.debug.(fn{ii}) = def.(fn{ii}); end, end
end

%% -------------------- LOAD -------------------- %%
tLOAD = dbg_tic('LOAD: rasters + ancillary + sinks', cfg);

S = struct();
if run.load_data
  fprintf('Loading data (means) ...\n');
  t=load('datasets_for_gmin/GHF_Hazzard_interp.mat','GHF_Hazzard_interp');  GHF_Hazzard_interp=t.GHF_Hazzard_interp; clear t
  t=load('datasets_for_gmin/GHF_Martos_interp.mat','GHF_Martos_interp');    GHF_Martos_interp=t.GHF_Martos_interp;   clear t
  t=load('datasets_for_gmin/GHF_Shen_interp.mat','GHF_Shen_interp');        GHF_Shen_interp=t.GHF_Shen_interp;       clear t
  t=load('datasets_for_gmin/GHF_Stal_interp.mat','GHF_Stal_interp');        GHF_Stal_interp=t.GHF_Stal_interp;       clear t
  t=load('datasets_for_gmin/GHF_Losing_interp.mat','GHF_Losing_interp');    GHF_Losing_interp=t.GHF_Losing_interp;   clear t
  t=load('datasets_for_gmin/GHF_An_interp.mat','GHF_An_interp');            GHF_An_interp=t.GHF_An_interp;           clear t
  t=load('datasets_for_gmin/GHF_FoxMaule_interp.mat','GHF_FoxMaule_interp');GHF_FoxMaule_interp=t.GHF_FoxMaule_interp;clear t

  % Ice thickness + grid
  t = load('datasets_for_gmin/coldex_icethk.mat','Xgrid','Ygrid','H');
  S.Xgrid = t.Xgrid; S.Ygrid = t.Ygrid; S.H = t.H; clear t
  S.X_km = S.Xgrid/1000; S.Y_km = S.Ygrid/1000;
    
  % Sinks (use cfg.sinks.path/varname consistently)
  sv = cfg.sinks.varname;  % e.g., 'sink_mask'
  ts = load(cfg.sinks.path, sv);
  assert(isfield(ts, sv), 'Variable "%s" not found in %s', sv, cfg.sinks.path);
  S.sink_mask = logical(ts.(sv)); clear ts
    
  % Specularity
  t = load('datasets_for_gmin/specularity.mat','Q'); 
  S.Q = t.Q; clear t
  if cfg.debug.quicklook_all, quicklook_field('Q (specularity)', S.X_km, S.Y_km, S.Q, S.sink_mask, cfg); end
    
  % Surface temp
  t = load('datasets_for_gmin/Ts_interp.mat','Ts_interp'); 
  S.Ts = t.Ts_interp; clear t
  if cfg.debug.quicklook_all, quicklook_field('Ts (surface T)', S.X_km, S.Y_km, S.Ts, [], cfg); end

  t=load('datasets_for_gmin/Mb_interp.mat','Mb_interp');                        S.Mb=t.Mb_interp; clear t
  %S.Mb = S.Mb * 0.6;
  if cfg.debug.quicklook_all, quicklook_field('Mb (accumulation)', S.X_km, S.Y_km, S.Mb, [], cfg); end

  t=load('datasets_for_gmin/mouginot_icevel.mat','icevel');          S.icevel=t.icevel; clear t
  if cfg.debug.quicklook_all, quicklook_field('Mb (accumulation)', S.X_km, S.Y_km, S.Mb, [], cfg); end

  % ---- Orientation sanity: ensure X and Y increase, flip ALL rasters coherently ----
  needFlipY = numel(S.Y_km)>=2 && S.Y_km(2,1) < S.Y_km(1,1);
  needFlipX = numel(S.X_km)>=2 && S.X_km(1,2) < S.X_km(1,1);    
  flipY = @(A) flipud(A);
  flipX = @(A) fliplr(A);
    
  baseList = {'H','Q','Ts','Mb','icevel','sink_mask'};
  for f = baseList
    S.(f{1}) = flip_if(S.(f{1}), needFlipY, needFlipX);
  end
  if needFlipY, S.Y_km = flipY(S.Y_km); end
  if needFlipX, S.X_km = flipX(S.X_km); end
    
  S.didFlipY = needFlipY;
  S.didFlipX = needFlipX;

  assert(isequal(size(S.H), size(S.Q), size(S.icevel), size(S.X_km), size(S.Y_km), size(S.sink_mask), size(S.Mb), size(S.Ts)), ...
    'All rasters must have identical sizes.');

  % Models (+ titles)
  S.models = struct('FoxMaule',GHF_FoxMaule_interp,'Hazzard',GHF_Hazzard_interp,'Martos',GHF_Martos_interp,'Shen',GHF_Shen_interp,'Stal',GHF_Stal_interp,'An',GHF_An_interp,'Losing',GHF_Losing_interp);
  S.names  = {'FoxMaule','Hazzard','Martos','Shen','Stal','An','Losing'};
  S.titles = {'Fox Maule','Hazzard','Martos','Shen','Stål','An','Lösing'};

  % Make sure model rasters match grid orientation
  if isfield(S,'didFlipY') && (S.didFlipY || S.didFlipX)
      fn = fieldnames(S.models);
      for k = 1:numel(fn)
          S.models.(fn{k}) = flip_if(S.models.(fn{k}), S.didFlipY, S.didFlipX);
      end
  end
 
  % --- Append flat GHF sweep models ---
  if isfield(cfg,'flat') && cfg.flat.enable
      fv = cfg.flat.values(:)';  % row
      for v = fv
          fld = matlab.lang.makeValidName(sprintf('Flat_%0.4g', v));
          S.models.(fld) = ones(size(S.H), 'like', S.H) * v;
          S.names{end+1}  = fld;
          S.titles{end+1} = sprintf('Flat %0.4g', v);
          if ~isfield(S,'unc'), S.unc = struct(); end
          S.unc.(fld) = zeros(size(S.H), 'like', S.H); % flat: σ=0
      end
  end

  % Uncertainty/bounds
  fprintf('Loading uncertainty/bounds ...\n');
  S.unc = struct(); S.bounds = struct();
  if isfile('datasets_for_gmin/UNC_Hazzard_interp.mat'), t=load('datasets_for_gmin/UNC_Hazzard_interp.mat','UNC_Hazzard_interp'); S.unc.Hazzard=t.UNC_Hazzard_interp; clear t; end
  if isfile('datasets_for_gmin/UNC_Martos_interp.mat'),   t=load('datasets_for_gmin/UNC_Martos_interp.mat','UNC_Martos_interp');   S.unc.Martos=t.UNC_Martos_interp;   clear t; end
  if isfile('datasets_for_gmin/UNC_Shen_interp.mat'),     t=load('datasets_for_gmin/UNC_Shen_interp.mat','UNC_Shen_interp');       S.unc.Shen=t.UNC_Shen_interp;       clear t; end
  if isfile('datasets_for_gmin/UNC_Stal_interp.mat'),     t=load('datasets_for_gmin/UNC_Stal_interp.mat','UNC_Stal_interp');       S.unc.Stal=t.UNC_Stal_interp;       clear t; end
  if isfile('datasets_for_gmin/BMIN_Losing_interp.mat'),  t=load('datasets_for_gmin/BMIN_Losing_interp.mat','BMIN_Losing_interp'); S.bounds.Losing_min=t.BMIN_Losing_interp; clear t; end
  if isfile('datasets_for_gmin/BMAX_Losing_interp.mat'),  t=load('datasets_for_gmin/BMAX_Losing_interp.mat','BMAX_Losing_interp'); S.bounds.Losing_max=t.BMAX_Losing_interp; clear t; end
    
  % Align uncertainties/bounds too
  if isfield(S,'didFlipY') && (S.didFlipY || S.didFlipX)
      if isfield(S,'unc') && ~isempty(S.unc)
          fu = fieldnames(S.unc);
          for k = 1:numel(fu)
              if ~isempty(S.unc.(fu{k})), S.unc.(fu{k}) = flip_if(S.unc.(fu{k}), S.didFlipY, S.didFlipX); end
          end
      end
      if isfield(S,'bounds')
          if isfield(S.bounds,'Losing_min') && ~isempty(S.bounds.Losing_min)
              S.bounds.Losing_min = flip_if(S.bounds.Losing_min, S.didFlipY, S.didFlipX);
          end
          if isfield(S.bounds,'Losing_max') && ~isempty(S.bounds.Losing_max)
              S.bounds.Losing_max = flip_if(S.bounds.Losing_max, S.didFlipY, S.didFlipX);
          end
      end
  end

  if cfg.debug.quicklook_all
      % Models
      fn = fieldnames(S.models);
      for k = 1:numel(fn)
          if startsWith(fn{k}, 'Flat_')
              continue;   % skip quicklook for flat models
          end
          quicklook_field(['Model: ' fn{k}], S.X_km, S.Y_km, S.models.(fn{k}), S.sink_mask, cfg);
      end

      % Uncertainty grids (if present)
      if isfield(S,'unc') && ~isempty(S.unc)
          fu = fieldnames(S.unc);
          for k = 1:numel(fu)
              if ~isempty(S.unc.(fu{k})), quicklook_field(['σ: ' fu{k}], S.X_km, S.Y_km, S.unc.(fu{k}), S.sink_mask, cfg); end
          end
      end

        % Losing bounds (if present)
        if isfield(S,'bounds')
            if isfield(S.bounds,'Losing_min') && ~isempty(S.bounds.Losing_min) 
                quicklook_field('Losing min', S.X_km, S.Y_km, S.bounds.Losing_min, S.sink_mask, cfg); end
            if isfield(S.bounds,'Losing_max') && ~isempty(S.bounds.Losing_max)
                quicklook_field('Losing max', S.X_km, S.Y_km, S.bounds.Losing_max, S.sink_mask, cfg); end
        end
  end  
end
%% -------------------- END LOAD -------------------- %%

% Cast → single (saves memory)
if cfg.to_single
  S.H = single(S.H); S.Q = single(S.Q); S.icevel = single(S.icevel);
  S.Ts = single(S.Ts); S.Mb = single(S.Mb);
  f = fieldnames(S.models); for k = 1:numel(f), S.models.(f{k}) = single(S.models.(f{k})); end
  if isfield(S,'unc'), f=fieldnames(S.unc); for k=1:numel(f), if ~isempty(S.unc.(f{k})), S.unc.(f{k}) = single(S.unc.(f{k})); end, end, end
  if isfield(S,'bounds')
  if isfield(S.bounds,'Losing_min') && ~isempty(S.bounds.Losing_min), S.bounds.Losing_min = single(S.bounds.Losing_min); end
  if isfield(S.bounds,'Losing_max') && ~isempty(S.bounds.Losing_max), S.bounds.Losing_max = single(S.bounds.Losing_max); end
  end
end

dbg_toc(tLOAD, 'LOAD: rasters + ancillary + sinks', cfg);

% Raster quick sanity
print_grid_stats(cfg, 'H (thickness)', S.H);
print_grid_stats(cfg, 'Q (specularity)', S.Q);
print_grid_stats(cfg, 'icevel', S.icevel);
print_grid_stats(cfg, 'sink_mask', S.sink_mask);

dbg_checkpoint('post_load', struct('X_km',S.X_km,'Y_km',S.Y_km,'H',S.H,'sink_mask',S.sink_mask), cfg);

%% -------------------- MATH CHECKS: post-load --------------------
assert_all_finite(S.X_km, 'X_km');  assert_all_finite(S.Y_km, 'Y_km');
assert_all_finite(S.H, 'H');        assert_all_finite(S.Q, 'Q');
negH = S.H(:) <= -1;
assert(~any(negH), 'H has values <= -1 (log(H+1) invalid)');
assert_in_range(S.Q, -10, 10, 'Q (specularity)');
badV = isfinite(S.icevel(:)) & (S.icevel(:) < 0);
assert(~any(badV), 'icevel has negative values');
assert(islogical(S.sink_mask), 'sink_mask must be logical');

%% -------------------- ROI MASK (sinks restrict ROI; with connectivity) ---------
fast_mask  = isfinite(S.icevel) & (S.icevel > cfg.v_keep);
west_mask  = (S.X_km >= 0) & (S.Y_km <= 350);

% ROI is everything slow+west+valid+in sinks
S.roi_mask = ~fast_mask & west_mask & isfinite(S.Q) & isfinite(S.H) & S.sink_mask;

idx        = find(S.roi_mask);
y_raw_vec  = S.Q(idx) >= cfg.spec_thresh;
H_vec      = S.H(idx);

% --- Connected components of ROI pixels ---
conn = cfg.sinks.connectivity;   % 4 or 8
CC   = bwconncomp(S.roi_mask, conn);

dbg(1,cfg,'ROI: %d pixels (%.1f%% of grid)', nnz(S.roi_mask), 100*nnz(S.roi_mask)/numel(S.roi_mask));
dbg(2,cfg,'Connected components (sinks) = %d', CC.NumObjects);

if cfg.debug.enable
    comp_sizes = cellfun(@numel, CC.PixelIdxList);
    figure('Name','Component-size histogram','Color','w');
    histogram(comp_sizes, 'BinWidth', 50);
    xlabel('pixels per component'); ylabel('count'); grid on;
    xlim([0 800]); ylim([-200 350]);
end

S.comp_id  = zeros(size(S.roi_mask),'uint32');
for k = 1:CC.NumObjects
    S.comp_id(CC.PixelIdxList{k}) = k;
end

try
    % labels inside ROI only (0 = non-ROI)
    Lc = S.comp_id;  Lc(~S.roi_mask) = 0;

    % ---- orientation-robust display copies ----
    Ydesc = numel(S.Y_km) > 1 && S.Y_km(2,1) < S.Y_km(1,1);
    if Ydesc
        Hplot  = flipud(S.H);
        Lcplot = flipud(Lc);
        Yplot  = flipud(S.Y_km);
    else
        Hplot  = S.H;
        Lcplot = Lc;
        Yplot  = S.Y_km;
    end
    
    Xplot = S.X_km;

    % ---- figure with two overlaid axes (separate colormaps) ----
    figC = figure('Color','w','Position',[100 100 950 820]);

    % background axes (H, grayscale)
    axH = axes('Parent',figC); hold(axH,'on'); axis(axH,'image');
    imagesc(axH, Xplot(1,:), Yplot(:,1), Hplot);
    set(axH,'YDir','normal','Layer','bottom'); box(axH,'on');
    colormap(axH, gray);
    try   % gentle contrast inside ROI
        Hroi = double(Hplot(S.roi_mask));
        if ~isempty(Hroi), clim(axH, prctile(Hroi,[2 98])); end
    catch
    end
    xlabel(axH,'X (km)'); ylabel(axH,'Y (km)');

    % overlay axes (components, categorical colormap)
    axC = axes('Parent',figC,'Color','none'); hold(axC,'on'); axis(axC,'image');
    set(axC,'YDir','normal','Layer','top'); box(axC,'on');
    hComp = imagesc(axC, Xplot(1,:), Yplot(:,1), Lcplot);
    set(hComp, 'AlphaData', 0.35 * double(Lcplot>0));
    % a colormap that includes index 0 (background) + components
    nC = double(max(Lcplot(:)));
    Cmap = lines(max(nC+1,2)); colormap(axC, Cmap);
    clim(axC, [0 max(nC,1)]);    % ensure indices map cleanly

    % ROI outline on top
    contour(axC, Xplot, Yplot, S.roi_mask, [0.5 0.5], 'k-', 'LineWidth', 0.8);

    % sync axes
    linkaxes([axH, axC]);
    xlim([0 800]); ylim([-200 350]);
    
    title(axC, sprintf('Connected sink components (n=%d)', nC));
    try exportgraphics(figC, fullfile(cfg.outdir,'map_components_over_H.png'), 'Resolution',300); catch, end
    if isfield(cfg,'debug') && isstruct(cfg.debug) && cfg.debug.enable && ...
       isfield(cfg.debug,'pause_on_fig') && cfg.debug.pause_on_fig
    end
catch
end

% ---- MATH CHECKS: ROI & components ----
assert(all(size(S.comp_id) == size(S.roi_mask)), 'comp_id size mismatch');
if CC.NumObjects > 0
    cvec = S.comp_id(S.roi_mask);
    assert(all(cvec(:) >= 1 & cvec(:) <= CC.NumObjects), 'comp_id out of 1..n_sinks over ROI');
end
assert(islogical(y_raw_vec) && numel(y_raw_vec)==nnz(S.roi_mask), 'y_raw_vec must be logical ROI-length');

dbg_checkpoint('post_roi', struct('roi_mask',S.roi_mask,'comp_id',S.comp_id), cfg);

comp_id_vec   = S.comp_id(idx);   % per-pixel component IDs
S.n_sinks     = CC.NumObjects;

% Cache dir keyed to ROI + sinks fingerprint (needed below)
roi_fp  = md5_of_logical(S.roi_mask);
sink_fp = md5_of_logical(S.sink_mask);
cfg.cache_dir = fullfile(cfg.outdir, sprintf('roi_cache_%s_%s', roi_fp(1:12), sink_fp(1:12)));
if ~exist(cfg.cache_dir,'dir'), mkdir(cfg.cache_dir); end

%% -------------------- LABELS (RAW + THK-ADJ) --------------------
if run.build_labels
  fprintf('Fitting RCS GLM and computing z-residuals...\n');
  Hz   = zscore(log(double(H_vec) + 1));
  q    = quantile(Hz, [0.05 0.275 0.50 0.725 0.95]);
  Xrcs = rcs_basis(Hz, q);
  [H_glm,~,~] = glmfit(Xrcs, [double(y_raw_vec) ones(size(y_raw_vec))], 'binomial', 'link','logit');
  p_thk     = glmval(H_glm, Xrcs, 'logit');
  z_resid   = (double(y_raw_vec) - p_thk) ./ sqrt(max(p_thk .* (1 - p_thk), eps));

  cands     = 0.00:0.05:3.50;
  pos_rate  = arrayfun(@(t) mean(z_resid > t), cands);
  [~,best]  = min(abs(100*pos_rate - cfg.target_pos_pct));
  z_thresh  = cands(best);
  fprintf('[Auto] z_thresh=%.2f (~%.1f%% adj positives)\n', z_thresh, 100*pos_rate(best));

  y_adj_vec = (z_resid > z_thresh);

  S.y_raw_vec = y_raw_vec; 
  S.y_adj_vec = y_adj_vec; 
  S.z_thresh  = z_thresh;
  S.y_adj_full = false(size(S.Q)); S.y_adj_full(idx) = y_adj_vec;
  S.y_raw_full = false(size(S.Q)); S.y_raw_full(idx) = y_raw_vec;

  fprintf('RAW positives: %d / %d (%.2f%%)\n', sum(y_raw_vec), numel(y_raw_vec), 100*mean(y_raw_vec));
  fprintf('ADJ positives: %d / %d (%.2f%%)\n', sum(y_adj_vec), numel(y_adj_vec), 100*mean(y_adj_vec));

  % ---- MATH CHECKS: labels ----
  assert_in_range(p_thk, 0, 1, 'p_thk (GLM prob)');
  assert_all_finite(z_resid, 'z_resid');
  adj_pct = 100*mean(y_adj_vec);
  if abs(adj_pct - cfg.target_pos_pct) > 2.5
      math_warn(cfg, 'Adjusted positive pct %.2f differs from target %.2f by >2.5%%', adj_pct, cfg.target_pos_pct);
  end
end

%% -------------------- PROCESS GHF → ROI cache --------------------
nM = numel(S.names);
gdif_cache = cell(nM,1);
S.model_cvals = nan(nM,1);

for i = 1:nM
    fprintf('[Uncertainty] Model %d/%d: %s\n', i, nM, S.names{i});
    dbg_progress(i, nM, cfg, cfg.debug.progress_every);
    name_i = S.names{i};
    Mi = S.models.(name_i);
    if startsWith(name_i,'Flat_')
        v = sscanf(name_i,'Flat_%f'); S.model_cvals(i) = v;
    else
        S.model_cvals(i) = median(Mi(idx),'omitnan');
    end
end

fprintf('Processing models to Gdif (ROI cache)...\n');

% Gmin (one-time)
gmin_path = fullfile(cfg.cache_dir, 'Gmin_full.mat');
if isfile(gmin_path)
    L = load(gmin_path, 'Gmin_full'); Gmin_full = L.Gmin_full; clear L
else
    [~, Gmin_full, ~] = processGHF(S.models.Hazzard, S.Ts, S.Mb, S.H); % model choice irrelevant for Gmin
    try save(gmin_path, 'Gmin_full', '-v7.3'); catch, warning('Could not cache Gmin_full'); end
end
gmin_vec = Gmin_full(idx);

if cfg.debug.quicklook_all
    quicklook_field('Gmin (baseline)', S.X_km, S.Y_km, Gmin_full, S.roi_mask, cfg);
end

% Per-model GΔ (ROI vectors) + cache
for i = 1:nM
   cache_i = fullfile(cfg.cache_dir, sprintf('gdifvec_%s.mat', S.names{i}));
   if isfile(cache_i)
       L = load(cache_i, 'gdif_vec'); gdif_cache{i} = L.gdif_vec;
   else
       Mi = S.models.(S.names{i});
       if startsWith(S.names{i}, 'Flat_') || (isscalar(unique(Mi(~isnan(Mi)))))
           gd = Mi(idx) - gmin_vec;
       else
          t1 = dbg_tic(sprintf('processGHF -> GΔ (%s)', S.names{i}), cfg);
          try
              [~, ~, tmpG] = processGHF(Mi, S.Ts, S.Mb, S.H);
          catch ME
              dbg(1,cfg,'ERROR processGHF (%s): %s', S.names{i}, ME.message);
              disp(getReport(ME,'extended')); rethrow(ME);
          end
          dbg_toc(t1, sprintf('processGHF -> GΔ (%s)', S.names{i}), cfg);
          gd = tmpG(idx);
       end
       gdif_cache{i} = gd; try gdif_vec = gd; save(cache_i, 'gdif_vec', '-v7.3'); catch, end
   end
end

% Build a single evaluation mask (apples-to-apples across models)
roiN = numel(idx);
valid_all = isfinite(gmin_vec(:));
for i = 1:nM
    vi = gdif_cache{i};
    if ~isvector(vi) || numel(vi) ~= roiN
        Mi = S.models.(S.names{i});
        if startsWith(S.names{i},'Flat_') || (isscalar(unique(Mi(~isnan(Mi)))))
            vi = Mi(idx) - gmin_vec;
        else
            [~, ~, Gdif_full] = processGHF(Mi, S.Ts, S.Mb, S.H);
            vi = Gdif_full(idx);
        end
        gdif_cache{i} = vi;
        try
            cache_i = fullfile(cfg.cache_dir, sprintf('gdifvec_%s.mat', S.names{i}));
            gdif_vec = vi(:); 
            save(cache_i, 'gdif_vec','-v7.3');
        catch
        end
    end
    vi = vi(:);
    valid_all = valid_all & isfinite(vi);
end

S.eval_mask = valid_all;

% Build 2D evaluation mask (full grid) from ROI-vector mask
S.eval_mask_full = false(size(S.roi_mask));
S.eval_mask_full(idx(S.eval_mask)) = true;

fprintf('Eval pixels: %d ; ROI pixels: %d ; ratio=%.3f\n', ...
    nnz(S.eval_mask), nnz(S.roi_mask), nnz(S.eval_mask)/nnz(S.roi_mask));

dbg(1,cfg,'Eval mask: %d of %d ROI pixels are valid across all models', nnz(S.eval_mask), numel(S.eval_mask));
dbg_checkpoint('post_eval_mask', struct('eval_mask',S.eval_mask), cfg);

% ---- MATH CHECKS: eval mask coherence ----
M_eval = S.eval_mask;
assert(nnz(M_eval)>0, 'eval_mask is empty');
assert(numel(M_eval)==numel(gmin_vec), 'eval_mask length mismatch');

for i=1:nM, assert(numel(gdif_cache{i})==numel(gmin_vec), 'gdif_cache len mismatch'); end

%% === Bootstrap groups (connected sinks) precompute (for 'component' mode) ===
S.boot_groups = struct('nC',0,'groups',{{}},'uniqC',[]);

if strcmpi(cfg.uncertainty.bootstrap_mode, 'component')
    Mglob = S.eval_mask;                    % boolean mask on ROI-vector space
    compM = comp_id_vec(Mglob);             % per-pixel component IDs within eval mask
    if ~isempty(compM)
        [uniqC,~,gidx] = unique(compM);
        groups = accumarray(gidx, (1:numel(compM))', [], @(v){v});  % cell of pixel-index lists
        S.boot_groups.uniqC  = uniqC;
        S.boot_groups.groups = groups;
        S.boot_groups.nC     = numel(groups);
    end

    if strcmpi(cfg.uncertainty.bootstrap_mode, 'component')
        dbg(1,cfg,'Bootstrap groups: %d components in eval-mask space', S.boot_groups.nC);
    end
end

%% -------------------- Estimate σ(Gmin) via finite differences -------------
if cfg.gmin_unc.enable
    tFD = dbg_tic('Finite-diff sigma(Gmin)', cfg);
        [S.sigma_gmin, S.sigma_gmin_scalar] = estimate_sigma_gmin_fd(S, cfg);
        dbg_toc(tFD, 'Finite-diff sigma(Gmin)', cfg);
        dbg(1,cfg,'sigma(Gmin) ROI RMS ≈ %.2f mW m^-2', S.sigma_gmin_scalar);
        dbg_checkpoint('post_sigma_gmin', struct('sigma_gmin_scalar',S.sigma_gmin_scalar), cfg);
    else
    S.sigma_gmin = [];
    S.sigma_gmin_scalar = NaN;
end

% ---- MATH CHECKS: sigma_gmin ----
if ~isempty(S.sigma_gmin)
    assert_all_finite(S.sigma_gmin, 'sigma_gmin');
    assert(all(S.sigma_gmin(:) >= 0 | isnan(S.sigma_gmin(:))), 'sigma_gmin must be >= 0');
    if S.sigma_gmin_scalar < 0
        math_warn(cfg,'sigma_gmin_scalar < 0 (%.3f) — impossible', S.sigma_gmin_scalar);
    end
end

%% -------------------- EVALUATE (point metrics, no uncertainty) ------------
if run.evaluate
  fprintf('Evaluating point metrics on ROI...\n');

  TP_r  = zeros(nM,1,'double'); FP_r  = TP_r; TN_r = TP_r; FN_r = TP_r;
  TP_a  = TP_r; FP_a  = TP_r;   TN_a  = TP_r; FN_a = TP_r;
  ACC_r = TP_r; ACC_a = TP_r;   PR_r  = TP_r; PR_a = TP_r; RC_r = TP_r; RC_a = TP_r; F1_r = TP_r; F1_a = TP_r;
  TPFP_r= TP_r; TNFN_r= TP_r;   TPFP_a= TP_r; TNFN_a= TP_r;

  for i = 1:nM
      dbg_progress(i, nM, cfg, cfg.debug.progress_every);
      fprintf('[Evaluate] Model %d/%d: %s\n', i, nM, S.names{i});
      x = gdif_cache{i};
      Mv   = S.eval_mask;

      prd  = x(Mv) > 0;
      yr   = S.y_raw_vec(Mv);
      ya   = S.y_adj_vec(Mv);

      m_raw = compute_metrics(prd, yr);
      m_adj = compute_metrics(prd, ya);

      % ---- MATH CHECKS: metrics identities ----
      assert(check_confusion_identities(m_raw), 'Metrics identities failed (raw) for %s', S.names{i});
      assert(check_confusion_identities(m_adj), 'Metrics identities failed (adj) for %s', S.names{i});

      TP_r(i)=m_raw.TP; FP_r(i)=m_raw.FP; TN_r(i)=m_raw.TN; FN_r(i)=m_raw.FN;
      TP_a(i)=m_adj.TP; FP_a(i)=m_adj.FP; TN_a(i)=m_adj.TN; FN_a(i)=m_adj.FN;
      ACC_r(i)=m_raw.ACC; ACC_a(i)=m_adj.ACC; PR_r(i)=m_raw.PREC; PR_a(i)=m_adj.PREC;
      RC_r(i)=m_raw.REC;  RC_a(i)=m_adj.REC;  F1_r(i)=m_raw.F1;   F1_a(i)=m_adj.F1;
      TPFP_r(i)=m_raw.TP_FP; TNFN_r(i)=m_raw.TN_FN; TPFP_a(i)=m_adj.TP_FP; TNFN_a(i)=m_adj.TN_FN;

      if mod(i, max(2,round(nM/10)))==0
        dbg(2,cfg,'[Eval] sample prec/rec (raw/adj) = %.3f/%.3f ; %.3f/%.3f', ...
            PR_r(i), RC_r(i), PR_a(i), RC_a(i));
      end
  end

  S.results_table = table( ...
    S.titles(:), nnz(idx)*ones(nM,1), ...
    TP_r, FP_r, TN_r, FN_r, TP_a, FP_a, TN_a, FN_a, ...
    ACC_r, ACC_a, PR_r, PR_a, RC_r, RC_a, F1_r, F1_a, ...
    TPFP_r, TNFN_r, TPFP_a, TNFN_a, ...
    'VariableNames', {'Model','N_roi', ...
    'TP_raw','FP_raw','TN_raw','FN_raw','TP_adj','FP_adj','TN_adj','FN_adj', ...
    'ACC_raw','ACC_adj','PREC_raw','PREC_adj','REC_raw','REC_adj','F1_raw','F1_adj', ...
    'TP_FP_raw','TN_FN_raw','TP_FP_adj','TN_FN_adj'});

    dbg_checkpoint('post_point_metrics', struct('results_table',S.results_table), cfg);
end

%% ---- OPTION A: expected metrics + pixel-bootstrap CIs (generalized) ----
if isfield(cfg,'uncertainty') && strcmpi(cfg.uncertainty.mode,'analytic')
  fprintf('Analytic expected metrics with pixel-bootstrap CIs...\n');

  nM   = numel(S.names);
  y_raw = S.y_raw_full(idx);
  y_adj = S.y_adj_full(idx);
  haveG = ~isempty(S.sigma_gmin);
  if haveG, sigG = S.sigma_gmin(idx); else, sigG = 0; end

  % we'll compute these metrics for both modes (raw/adj)
  stat_list = {'SPEC','REC','PREC','ACC','F1'};

  % containers for "expected" (non-bootstrap) per-model metrics
  EXP_raw = struct(); EXP_adj = struct();
  for s = stat_list, EXP_raw.(s{1}) = nan(nM,1); EXP_adj.(s{1}) = nan(nM,1); end

  % containers for bootstrap CIs (outer/inner) per metric (nM x 2)
  CI_raw  = struct();  CI_adj  = struct();
  CIi_raw = struct();  CIi_adj = struct();
  for s = stat_list
      CI_raw.(s{1})  = nan(nM,2);  CI_adj.(s{1})  = nan(nM,2);
      CIi_raw.(s{1}) = nan(nM,2);  CIi_adj.(s{1}) = nan(nM,2);
  end

  alpha_main = 1 - cfg.uncertainty.band_main; % e.g., 0.10 for 90% 
  alpha_inner = 1 - cfg.uncertainty.band_inner; % e.g., 0.50 for 50%

  for i = 1:nM
      dbg_progress(i, nM, cfg, cfg.debug.progress_every);
      name_i = S.names{i};
      x  = gdif_cache{i};
      M  = S.eval_mask;
      if ~any(M), continue; end

      % per-pixel σ for model on ROI
      sigM = zeros(size(idx),'like',x);
      switch name_i
        case {'Hazzard','Martos','Shen','Stal'}
          if isfield(S.unc,name_i) && ~isempty(S.unc.(name_i))
              sigM = S.unc.(name_i)(idx);
          else
              sigM(:) = get_field_default(fallback_err, name_i, 0);
          end
        case 'Losing'
          if isfield(S.bounds,'Losing_min') && isfield(S.bounds,'Losing_max') ...
             && ~isempty(S.bounds.Losing_min) && ~isempty(S.bounds.Losing_max)
              half = 0.5*abs(S.bounds.Losing_max(idx) - S.bounds.Losing_min(idx));
              sigM = half / 1.96;    % approx σ from 95% range
          else
              sigM(:) = get_field_default(fallback_err, name_i, 0);
          end
        otherwise
          sigM(:) = get_field_default(fallback_err, name_i, 0);
      end

      % Pr[wet] = Pr(GΔ > 0) with combined σ (model ⊕ Gmin)
      pw = wet_prob_from_sigmas(x, sigM, haveG*sigG);

        % tighten to finite-only mask for this model
        I = M & isfinite(pw) & isfinite(y_raw) & isfinite(y_adj);
        
        if ~all(I(:) == M(:))
            dbg(1,cfg,'[%s] dropping %d/%d non-finite pixels from eval mask', ...
                name_i, sum(M(:) & ~I(:)), sum(M(:)));
        end
        
        % Use only finite entries from here on
        pwM = pw(I);
        yrM = double(y_raw(I));   % cast to double to avoid any logical*double quirks
        yaM = double(y_adj(I));
        
        % (optional but robust) clamp tiny numeric drift
        pwM = min(max(pwM,0),1);


      % ---- MATH CHECKS: probability + expectations ----
      assert_in_range(pw, 0, 1, sprintf('Pr(GΔ>0) pw (%s)', name_i));
      assert_all_finite(pw, 'pw');
      if any(isfinite(pwM))
          sigM2 = sigM * 1.5;
          pw2   = wet_prob_from_sigmas(x, sigM2, haveG*sigG);
          d = abs(pw2(M) - 0.5) - abs(pwM - 0.5);
          if mean(d,'omitnan') > 1e-3
              math_warn(cfg, '%s: pw did not move toward 0.5 with larger σ (mean Δ=%.3g)', name_i, mean(d,'omitnan'));
          end
      end

      % ----------------------- Expected (non-bootstrap) confusion and metrics ----------
      % RAW expected counts
      ETP  = sum(pwM .* yrM);
      EFP  = sum(pwM .* (1-yrM));
      ETN  = sum((1-pwM) .* (1-yrM));
      EFN  = sum((1-pwM) .* yrM);
      mR   = metrics_from_counts(ETP,EFP,ETN,EFN);
      EXP_raw.SPEC(i) = mR.SPEC; EXP_raw.REC(i) = mR.REC; EXP_raw.PREC(i)=mR.PREC;
      EXP_raw.ACC(i)  = mR.ACC;  EXP_raw.F1(i)  = mR.F1;

      % ADJ expected counts
      ETPa = sum(pwM .* yaM);
      EFPa = sum(pwM .* (1-yaM));
      ETNa = sum((1-pwM) .* (1-yaM));
      EFNa = sum((1-pwM) .* yaM);
      mA   = metrics_from_counts(ETPa,EFPa,ETNa,EFNa);
      EXP_adj.SPEC(i) = mA.SPEC; EXP_adj.REC(i) = mA.REC; EXP_adj.PREC(i)=mA.PREC;
      EXP_adj.ACC(i)  = mA.ACC;  EXP_adj.F1(i)  = mA.F1;

      % skip bootstrapping for flat models
      if startsWith(name_i,'Flat_')
          % CI arrays were preinitialized to NaN, so just skip.
          continue
      end

      % ----------------------- Pixel bootstrap ---------------------------
      N = numel(pwM);
      if N==0
          for s = stat_list
              CI_raw.(s{1})(i,:)  = [NaN NaN]; CI_adj.(s{1})(i,:)  = [NaN NaN];
              CIi_raw.(s{1})(i,:) = [NaN NaN]; CIi_adj.(s{1})(i,:) = [NaN NaN];
          end
      else
      % ----------------------- Bootstrap (pixel or component) -----------------------
      B = cfg.uncertainty.n_boot;
      rng(42);
        
      % Allocate one vector per metric per mode (RAW/ADJ)
      BR = struct(); BA = struct();
      for s = stat_list, BR.(s{1}) = zeros(B,1); BA.(s{1}) = zeros(B,1); end
    
      modeStr = char(cfg.uncertainty.bootstrap_mode);
    
      if strcmpi(modeStr, 'pixel')
          % --- Per-pixel bootstrap (classic) ---
          for b = 1:B
              rb = randi(N, N, 1);   % resample pixels with replacement
              % RAW expected counts
              ETPb = sum(pwM(rb) .* yrM(rb));
              EFPb = sum(pwM(rb) .* (1 - yrM(rb)));
              ETNb = sum((1 - pwM(rb)) .* (1 - yrM(rb)));
              EFNb = sum((1 - pwM(rb)) .* yrM(rb));
              mbR  = metrics_from_counts(ETPb,EFPb,ETNb,EFNb);
              % ADJ expected counts
              ETPab = sum(pwM(rb) .* yaM(rb));
              EFPab = sum(pwM(rb) .* (1 - yaM(rb)));
              ETNab = sum((1 - pwM(rb)) .* (1 - yaM(rb)));
              EFNab = sum((1 - pwM(rb)) .* yaM(rb));
              mbA   = metrics_from_counts(ETPab,EFPab,ETNab,EFNab);
              % collect
              BR.SPEC(b)=mbR.SPEC; BR.REC(b)=mbR.REC; BR.PREC(b)=mbR.PREC; BR.ACC(b)=mbR.ACC; BR.F1(b)=mbR.F1;
              BA.SPEC(b)=mbA.SPEC; BA.REC(b)=mbA.REC; BA.PREC(b)=mbA.PREC; BA.ACC(b)=mbA.ACC; BA.F1(b)=mbA.F1;
          end
                        
              if cfg.debug.enable
                  fprintf('[CI diag, pixel boot] %-10s  Npix=%d  nBoot=%d  mu(pw)=%.3f  sd(pw)=%.3f  frac~{0,1}=%.2f\n', ...
                          name_i, N, B, mean(pwM,'omitnan'), std(pwM,'omitnan'), ...
                          mean((pwM<1e-3)|(pwM>1-1e-3),'omitnan'));
              end        
          
      elseif strcmpi(modeStr, 'component')
          % --- Sink-component bootstrap (resample connected components) ---
          if ~isfield(S,'boot_groups') || S.boot_groups.nC == 0
              % Fallback if no groups available
              for s = stat_list
                  CI_raw.(s{1})(i,:)  = [NaN NaN]; CI_adj.(s{1})(i,:)  = [NaN NaN];
                  CIi_raw.(s{1})(i,:) = [NaN NaN]; CIi_adj.(s{1})(i,:) = [NaN NaN];
              end
              % Skip quantiles for this model
              continue
          end
        
          G  = S.boot_groups.groups;     % cell of index vectors into M-space
          nC = S.boot_groups.nC;
    
          for b = 1:B
                rbC = randi(nC, nC, 1);      % sample components with replacement
                pick = vertcat(G{rbC});      % concatenate member pixels
    
                % RAW expected counts
                ETPb = sum(pwM(pick) .* yrM(pick));
                EFPb = sum(pwM(pick) .* (1 - yrM(pick)));
                ETNb = sum((1 - pwM(pick)) .* (1 - yrM(pick)));
                EFNb = sum((1 - pwM(pick)) .* yrM(pick));
                mbR  = metrics_from_counts(ETPb,EFPb,ETNb,EFNb);
    
                % ADJ expected counts
                ETPab = sum(pwM(pick) .* yaM(pick));
                EFPab = sum(pwM(pick) .* (1 - yaM(pick)));
                ETNab = sum((1 - pwM(pick)) .* (1 - yaM(pick)));
                EFNab = sum((1 - pwM(pick)) .* yaM(pick));
                mbA   = metrics_from_counts(ETPab,EFPab,ETNab,EFNab);
    
                % collect
                BR.SPEC(b)=mbR.SPEC; BR.REC(b)=mbR.REC; BR.PREC(b)=mbR.PREC; BR.ACC(b)=mbR.ACC; BR.F1(b)=mbR.F1;
                BA.SPEC(b)=mbA.SPEC; BA.REC(b)=mbA.REC; BA.PREC(b)=mbA.PREC; BA.ACC(b)=mbA.ACC; BA.F1(b)=mbA.F1;
           end            
               
           if cfg.debug.enable
               fprintf('[CI diag, comp boot]   %-10s  nComp=%d  nBoot=%d  mu(pw)=%.3f  sd(pw)=%.3f  frac~{0,1}=%.2f\n', ...
               name_i, nC, B, mean(pwM,'omitnan'), std(pwM,'omitnan'), ...
               mean((pwM<1e-3)|(pwM>1-1e-3),'omitnan'));
           end
       else
           error('cfg.uncertainty.bootstrap_mode must be ''pixel'' or ''component''.');
       end
        
       % Quantiles → outer/inner CI bands
       for s = stat_list
           v = BR.(s{1});
           CI_raw.(s{1})(i,:)  = quantile(v, [alpha_main/2, 1 - alpha_main/2]);
           CIi_raw.(s{1})(i,:) = quantile(v, [alpha_inner/2, 1 - alpha_inner/2]);
           v = BA.(s{1});
           CI_adj.(s{1})(i,:)  = quantile(v, [alpha_main/2, 1 - alpha_main/2]);
           CIi_adj.(s{1})(i,:) = quantile(v, [alpha_inner/2, 1 - alpha_inner/2]);

           % ---- MATH CHECKS: CI ordering & nesting ----
           lohiR  = CI_raw.(s{1})(i,:);  lohiRi = CIi_raw.(s{1})(i,:);
           lohiA  = CI_adj.(s{1})(i,:);  lohiAi = CIi_adj.(s{1})(i,:);
           assert(lohiR(1) <= lohiR(2)+1e-12,  '%s RAW outer CI inverted', s{1});
           assert(lohiRi(1)<= lohiRi(2)+1e-12, '%s RAW inner CI inverted', s{1});
           assert(lohiA(1) <= lohiA(2)+1e-12,  '%s ADJ outer CI inverted', s{1});
           assert(lohiAi(1)<= lohiAi(2)+1e-12, '%s ADJ inner CI inverted', s{1});
           if all(isfinite(lohiR)) && all(isfinite(lohiRi))
               if ~(lohiR(1)-1e-6 <= lohiRi(1) && lohiRi(2) <= lohiR(2)+1e-6)
                   math_warn(cfg,'%s RAW inner CI not within outer CI (ok in rare cases)', s{1});
               end
           end
           if all(isfinite(lohiA)) && all(isfinite(lohiAi))
               if ~(lohiA(1)-1e-6 <= lohiAi(1) && lohiAi(2) <= lohiA(2)+1e-6)
                   math_warn(cfg,'%s ADJ inner CI not within outer CI (ok in rare cases)', s{1});
               end
           end
        end
      end
  end
end

% ---- Stash everything into S for plotting ----
S.Gstats = struct();

% RAW
for s = stat_list
 S.Gstats.([s{1} '_raw']) = struct( ...
   'EXP', EXP_raw.(s{1}), ...
   'CI',  CI_raw.(s{1}), ...
   'CIi', CIi_raw.(s{1}) );
end
% ADJ
for s = stat_list
 S.Gstats.([s{1} '_adj']) = struct( ...
   'EXP', EXP_adj.(s{1}), ...
   'CI',  CI_adj.(s{1}), ...
   'CIi', CIi_adj.(s{1}) );
end

% quick CSV summary for ADJ (edit as needed)
try
 Tsum = table( string(S.names(:)), string(S.titles(:)), ...
     S.Gstats.SPEC_adj.EXP(:), S.Gstats.SPEC_adj.CI(:,1), S.Gstats.SPEC_adj.CI(:,2), ...
     S.Gstats.REC_adj.EXP(:),  S.Gstats.REC_adj.CI(:,1),  S.Gstats.REC_adj.CI(:,2), ...
     S.Gstats.PREC_adj.EXP(:), S.Gstats.PREC_adj.CI(:,1), S.Gstats.PREC_adj.CI(:,2), ...
     S.Gstats.ACC_adj.EXP(:),  S.Gstats.ACC_adj.CI(:,1),  S.Gstats.ACC_adj.CI(:,2), ...
     S.Gstats.F1_adj.EXP(:),   S.Gstats.F1_adj.CI(:,1),   S.Gstats.F1_adj.CI(:,2), ...
     'VariableNames', {'Name','Title', ...
       'SPEC','SPEC_lo','SPEC_hi','REC','REC_lo','REC_hi', ...
       'PREC','PREC_lo','PREC_hi','ACC','ACC_lo','ACC_hi','F1','F1_lo','F1_hi'});
 if ~exist(cfg.outdir,'dir'), mkdir(cfg.outdir); end
 csv_path = fullfile(cfg.outdir,'ci_summary_adj.csv');
 writetable(Tsum, csv_path);
 fprintf('Saved CI summary CSV (ADJ): %s\n', csv_path);
catch ME
 warning(ME.identifier,'Could not write ci_summary_adj.csv: %s', ME.message);
end

dbg_checkpoint('post_bootstrap', struct('Gstats',S.Gstats), cfg);

%% -------------------- COLORED-CI PLOTS (single & generalized) ------------
if isfield(cfg,'two_stat_pairs') && ~isempty(cfg.two_stat_pairs)
    % ---- MATH CHECKS: plotting prerequisites ----
    assert(isfield(S,'Gstats'), 'Gstats missing — run uncertainty/bootstraps first');
    req = {'SPEC_adj','REC_adj','PREC_adj','ACC_adj','F1_adj'};
    have = isfield(S.Gstats, req);
    if ~all(have)
        miss = req(~have);
        math_warn(cfg, 'Missing Gstats fields for ADJ plots: %s', strjoin(miss, ', '));
    end

    for k = 1:numel(cfg.two_stat_pairs)
        pair = cfg.two_stat_pairs{k};
        try
            % Basic field presence check for this pair
            px = upper(pair{1}); py = upper(pair{2});
            needX = [px '_adj']; needY = [py '_adj'];
            if ~isfield(S.Gstats, needX) || ~isfield(S.Gstats, needY)
                math_warn(cfg, 'Skipping plot %s-%s (ADJ): missing %s or %s', px, py, needX, needY);
                continue;
            end

            % Extra sanity: CI order & finite EXP (vectorized quick check)
            GX = S.Gstats.(needX); GY = S.Gstats.(needY);
            if isfield(GX,'CI')
                badX = any(GX.CI(:,1) > GX.CI(:,2));
                if any(badX), math_warn(cfg,'%s ADJ outer CI inverted for some models', px); end
            end
            if isfield(GY,'CI')
                badY = any(GY.CI(:,1) > GY.CI(:,2));
                if any(badY), math_warn(cfg,'%s ADJ outer CI inverted for some models', py); end
            end
            assert(all(isfinite(GX.EXP) | isnan(GX.EXP)), '%s ADJ EXP has non-finite non-NaN entries', px);
            assert(all(isfinite(GY.EXP) | isnan(GY.EXP)), '%s ADJ EXP has non-finite non-NaN entries', py);

            % ---- Plot panel (ADJ) ----
            plot_colored_ci_pair(S, cfg, pair{1}, pair{2}, 'adj');
            out = fullfile(cfg.outdir, sprintf('ci2_%s_%s_ADJ.png', upper(pair{1}), upper(pair{2})));
            exportgraphics(gcf, out, 'Resolution', 300);
            if cfg.close_figs, close(gcf); end

            % If this pair is SPEC–REC, also save the legacy filename:
            if strcmpi(pair{1},'SPEC') && strcmpi(pair{2},'REC')
                exportgraphics(gcf, fullfile(cfg.outdir,'specrec_expected_ADJ.png'), 'Resolution',300);
            end
        catch ME
            warning('2D plot %s-%s ADJ failed: %s', pair{1}, pair{2}, ME.message);
        end
    end
end

% RAW plots for the same pairs
if isfield(cfg.plots,'prrec_expected_raw') && cfg.plots.prrec_expected_raw ...
        && isfield(cfg,'two_stat_pairs') && ~isempty(cfg.two_stat_pairs)
    for k = 1:numel(cfg.two_stat_pairs)
        pair = cfg.two_stat_pairs{k};
        try
            plot_colored_ci_pair(S, cfg, pair{1}, pair{2}, 'raw');
            out = fullfile(cfg.outdir, sprintf('ci2_%s_%s_RAW.png', upper(pair{1}), upper(pair{2})));
            exportgraphics(gcf, out, 'Resolution', 300);
            if cfg.close_figs, close(gcf); end
        catch ME
            warning('2D plot %s-%s RAW failed: %s', pair{1}, pair{2}, ME.message);
        end
    end
end

%% ===== DEBUG READOUTS (CI sanity + quick table) ===========================
if isfield(cfg,'debug') && isstruct(cfg.debug) && cfg.debug.enable
    try
        is_flat = startsWith(S.names,'Flat_')';
        not_flat = ~is_flat;
        Npix = nnz(S.eval_mask);

        % CI widths
        W_SPEC_raw = S.Gstats.SPEC_raw.CI(:,2) - S.Gstats.SPEC_raw.CI(:,1);
        W_REC_raw  = S.Gstats.REC_raw.CI(:,2)  - S.Gstats.REC_raw.CI(:,1);
        W_SPEC_adj = S.Gstats.SPEC_adj.CI(:,2) - S.Gstats.SPEC_adj.CI(:,1);
        W_REC_adj  = S.Gstats.REC_adj.CI(:,2)  - S.Gstats.REC_adj.CI(:,1);

        pr = @(v) [min(v,[],'omitnan') median(v,'omitnan') mean(v,'omitnan') max(v,[],'omitnan')];

        fprintf('\n========== CI / Bootstrap Diagnostics ==========\n');
        fprintf('Eval pixels (N): %d\n', Npix);
        fprintf('Models: %d total | non-flat: %d | flat: %d\n', numel(S.names), nnz(not_flat), nnz(is_flat));

        if any(not_flat)
            s1 = pr(W_SPEC_adj(not_flat)); s2 = pr(W_REC_adj(not_flat));
            r1 = pr(W_SPEC_raw(not_flat)); r2 = pr(W_REC_raw(not_flat));
            fprintf('\n[Adj] SPEC width  (non-flat)  min/med/mean/max = %.4f / %.4f / %.4f / %.4f\n', s1);
            fprintf('[Adj] REC  width  (non-flat)  min/med/mean/max = %.4f / %.4f / %.4f / %.4f\n', s2);
            fprintf('[Raw] SPEC width  (non-flat)  min/med/mean/max = %.4f / %.4f / %.4f / %.4f\n', r1);
            fprintf('[Raw] REC  width  (non-flat)  min/med/mean/max = %.4f / %.4f / %.4f / %.4f\n', r2);
        end

        % CSV summary (complete)
        Tsum = table( string(S.names(:)), string(S.titles(:)), ~is_flat(:), ...
            S.Gstats.SPEC_raw.EXP(:), S.Gstats.SPEC_raw.CI(:,1), S.Gstats.SPEC_raw.CI(:,2), W_SPEC_raw(:), ...
            S.Gstats.REC_raw.EXP(:),  S.Gstats.REC_raw.CI(:,1),  S.Gstats.REC_raw.CI(:,2),  W_REC_raw(:), ...
            S.Gstats.SPEC_adj.EXP(:), S.Gstats.SPEC_adj.CI(:,1), S.Gstats.SPEC_adj.CI(:,2), W_SPEC_adj(:), ...
            S.Gstats.REC_adj.EXP(:),  S.Gstats.REC_adj.CI(:,1),  S.Gstats.REC_adj.CI(:,2),  W_REC_adj(:), ...
            'VariableNames', {'Name','Title','IsModel', ...
                'SPEC_raw','SPEC_raw_lo','SPEC_raw_hi','W_SPEC_raw', ...
                'REC_raw','REC_raw_lo','REC_raw_hi','W_REC_raw', ...
                'SPEC_adj','SPEC_adj_lo','SPEC_adj_hi','W_SPEC_adj', ...
                'REC_adj','REC_adj_lo','REC_adj_hi','W_REC_adj'});

        if ~exist(cfg.outdir,'dir'), mkdir(cfg.outdir); end
        csv_path = fullfile(cfg.outdir,'ci_summary.csv');
        writetable(Tsum, csv_path);
        fprintf('\nSaved CI summary CSV: %s\n', csv_path);

    catch ME
        warning(ME.identifier,'CI/diagnostic readout failed: %s', ME.message);
    end
end

%% -------------------- MAPS -------------------------------------
if cfg.plots.maps
  fprintf('Rendering maps (GΔ masked by slow flow; sinks overlaid as points)...\n');
  tMAP = dbg_tic('Render maps', cfg);

  % ---- Masks & styles ----
  fast_mask = isfinite(S.icevel) & (S.icevel > cfg.v_keep);   % hide these
  slow_mask = isfinite(S.icevel) & ~fast_mask;                % show these

  sinkStyle = cfg.sinks;
  if ~isfield(sinkStyle,'marker'), sinkStyle.marker = '.'; end
  if ~isfield(sinkStyle,'size'),   sinkStyle.size   = 5;  end
  if ~isfield(sinkStyle,'color'),  sinkStyle.color  = [0.05 0.05 0.05]; end
  if ~isfield(sinkStyle,'alpha'),  sinkStyle.alpha  = 0.85; end
  sinkStyle.mode = 'points';   % force points, not fill

  % ---- Reference zero-centering for diverging clim (2–98%% over ROI) ----
  try
      [~, ~, Gdif_ref] = processGHF(S.models.Hazzard, S.Ts, S.Mb, S.H);
      Zref = double(Gdif_ref(S.roi_mask & isfinite(Gdif_ref)));
      if isempty(Zref), error('empty Zref'); end
      L = max(abs(prctile(Zref,[2 98])));
      if ~isfinite(L) || L <= 0, L = 50; end   % safe fallback
  catch
      warning('Could not compute reference L from Hazzard; using fallback L=50');
      L = 50;
  end
  fprintf('Diverging limits set to [-%.3g, %.3g]\n', L, L);

  % ---- GΔ maps for each non-flat model ----
  nM = numel(S.names);
  for ii = 1:nM
      name = S.names{ii};
      title_i = S.titles{ii};
      if startsWith(name,'Flat_')
          dbg(2,cfg,'[Maps] Skipping flat model: %s', name);
          continue;
      end

      % cache or compute full GΔ
      zd_path = fullfile(cfg.cache_dir, sprintf('Gdif_full_%s.mat', name));
      if isfile(zd_path)
          Ld = load(zd_path,'Gdif_full');  Zd = Ld.Gdif_full;
      else
          [~, ~, Zd] = processGHF(S.models.(name), S.Ts, S.Mb, S.H);
          try Gdif_full = Zd; save(zd_path,'Gdif_full','-v7.3'); catch, end
      end

      % plot: slow-flow alpha mask; sinks as points
      fig = plot_overlay_on_H( ...
          S.X_km, S.Y_km, S.H, Zd, [-L L], 'diverging', ...
          cfg, slow_mask, ...
          'SinkMask',  S.roi_mask, ...
          'SinkStyle', sinkStyle);

      title(sprintf('G_Δ: %s', title_i), 'FontSize', 14);
      xlabel('X (km)'); ylabel('Y (km)');
      xlim([0 800]);
      ylim([-200 350]);

      safeTitle = sanitize_filename(title_i);
      outfile   = fullfile(cfg.outdir, sprintf('map_Gdif_%s_over_H.png', safeTitle));
      ok = save_png(fig, outfile, cfg);
      if ~ok, warning('Map save failed for %s', title_i); end
      clear Zd Ld
  end

  % ---- G_min map (sequential) with same overlays ----
  try
      [~, Gmin_map, ~] = processGHF(S.models.Hazzard, S.Ts, S.Mb, S.H);
      [loGmin, hiGmin] = robust_clim_fast(Gmin_map, S.roi_mask);   % 2–98% over ROI

      fig = plot_overlay_on_H( ...
          S.X_km, S.Y_km, S.H, Gmin_map, [loGmin hiGmin], 'sequential', ...
          cfg, slow_mask, ...
          'SinkMask',  S.roi_mask, ...
          'SinkStyle', sinkStyle);

      title('G_{min} (baseline)','FontSize',14);
      xlabel('X (km)'); ylabel('Y (km)');
      xlim([0 800]);
      ylim([-200 350]);

      ok = save_png(fig, fullfile(cfg.outdir,'map_Gmin_over_H.png'), cfg);
      if ~ok, warning('Map save failed for Gmin'); end
  catch 
  end

  % ---- Final sanity ----
  assert(isfinite(L) && L>0, 'map bound L must be positive finite');

  dbg_toc(tMAP, 'Render maps', cfg);
  dbg_checkpoint('post_maps', struct('L',L), cfg);
end
% ===================== END MAPS =====================
if cfg.profile, profile viewer; profile off; end
end
% ===================== end main =====================

%% ===================== Helpers =====================
function A = flip_if(A, doY, doX)
% flip array A in Y and/or X if flags are true
    if doY, A = flipud(A); end
    if doX, A = fliplr(A); end
end

function dbg(level, cfg, fmt, varargin)
  if isfield(cfg,'debug') && cfg.debug.enable && cfg.debug.level>=level
      fprintf([fmt '\n'], varargin{:});
  end
end

function s = sanitize_filename(s)
% Replace anything non [A-Za-z0-9_.-] with _
    s = regexprep(char(s), '[^A-Za-z0-9_.-]+', '_');
end

function print_grid_stats(cfg, name, A)
  if ~cfg.debug.enable, return; end
  sz = size(A); f = isfinite(A); n=numel(A);
  p = @(v) prctile(double(v),[2 50 98]);
  if any(f(:))
      q = p(A(f));
      fprintf('[%s] size=%dx%d finite=%d/%.0f%% nan=%d min/med/max≈[%.3g %.3g %.3g]\n', name, sz(1), sz(2), nnz(f), 100*nnz(f)/n, n-nnz(f), q(1), q(2), q(3));
  else
      fprintf('[%s] size=%dx%d finite=0/%.0f%%\n', name, sz(1), sz(2), 100);
  end
end

function ok = save_png(fig, path, cfg)
% Robust PNG saver:
% - absolute path
% - settle graphics queue
% - try exportgraphics
% - if it fails (e.g., invalid/deleted rulers), clone to a headless fig and retry
% - final fallback to print

    if nargin<3 || isempty(cfg), cfg = struct('overwrite',true,'close_figs',false); end
    if ~ishandle(fig) || ~isvalid(fig)
        warning('save_png:badFig','Invalid figure handle'); ok=false; return;
    end

    % Ensure directory exists; force absolute path
    [pdir, pname, pext] = fileparts(path);
    if isempty(pext), pext = '.png'; end
    if ~exist(pdir,'dir'), mkdir(pdir); end
    apath = fullfile(char(java.io.File(pdir).getCanonicalPath), [pname pext]);

    if isfield(cfg,'overwrite') && ~cfg.overwrite && isfile(apath)
        fprintf('save_png: not overwriting existing %s\n', apath);
        ok = true;
        if isfield(cfg,'close_figs') && cfg.close_figs && isvalid(fig), close(fig); end
        return;
    end

    % Let any pending layout/colormap/CB updates finish
    drawnow expose;
    pause(0.01);

    ok = try_export(fig, apath);

    % If exportgraphics failed (often due to ruler/listener cleanup), retry on a clone
    if ~ok
        try
            fclone = clone_figure_for_export(fig);
            drawnow expose;
            pause(0.01);
            ok = try_export(fclone, apath);
            if isvalid(fclone), close(fclone); end
        catch MEc
            warning('save_png:cloneExport','Clone-export path failed: %s', MEc.message);
        end
    end

    % Final fallback: print
    if ~ok
        try
            if isvalid(fig), set(fig,'PaperPositionMode','auto'); end
            print(fig, apath, '-dpng','-r300');
            ok = true;
        catch ME2
            warning('save_png:print','print fallback failed: %s', ME2.message);
        end
    end

    if ok
        fprintf('Saved PNG → %s\n', apath);
    else
        warning('save_png:failed','Could not save PNG to %s', apath);
    end

    if isfield(cfg,'close_figs') && cfg.close_figs && isvalid(fig), close(fig); end
end

% ---- helpers for robust export ----
function ok = try_export(h, apath)
    ok = false;
    if ~ishandle(h) || ~isvalid(h), return; end
    try
        exportgraphics(h, apath, 'Resolution', 300);
        ok = true;
    catch ME
        warning('save_png:exportgraphics','exportgraphics failed: %s', ME.message);
    end
end

function f2 = clone_figure_for_export(f1)
% Create a headless clone of a figure to decouple export from live listeners/linkprops
    if ~ishandle(f1) || ~isvalid(f1)
        error('clone_figure_for_export:invalid','Source figure invalid.');
    end
    % Build new invisible figure roughly matching size/background
    pos = get(f1,'Position');
    try bg = get(f1,'Color'); catch, bg = [1 1 1]; end
    f2  = figure('Visible','off','Color',bg,'Position',pos,'InvertHardCopy','off');
    % Copy ALL children (axes, colorbars, legends, etc.) preserving stacking order
    kids = allchild(f1);
    copyobj(kids, f2);
    % Try to preserve axis limits & colormaps on topmost axes
    axs = findall(f2,'Type','axes');
    for k = 1:numel(axs)
        try axis(axs(k),'image'); catch, end
    end
end

function m = compute_metrics(pred, truth)
  pred  = logical(pred(:)); truth = logical(truth(:));
  N  = numel(truth);
  TP = sum( pred &  truth);
  FP = sum( pred & ~truth);
  TN = sum(~pred & ~truth);
  FN = sum(~pred &  truth);
  ACC  = (TP + TN) / max(N, eps);
  PREC = TP / max(TP + FP, eps);
  SPEC = TN / max(TN + FP, eps);
  REC  = TP / max(TP + FN, eps);
  F1   = 2 * PREC * REC / max(PREC + REC, eps);
  m = struct('N',N,'TP',TP,'FP',FP,'TN',TN,'FN',FN, ...
             'ACC',ACC,'PREC',PREC,'SPEC',SPEC,'REC',REC,'F1',F1, ...
             'TP_FP', TP / max(FP, eps), ...
             'TN_FN', TN / max(FN, eps));
end

function [xExt, yExt] = grid_extents(X, Y)
% Return [xmin xmax] and [ymin ymax] at pixel edges (robust to NaNs, order)
  % X spacing
  x = X(1,:); x = x(isfinite(x));
  if numel(x) >= 2
      dxv = diff(x); dxv = dxv(isfinite(dxv) & dxv~=0);
      if isempty(dxv), dx = 1; else, dx = median(dxv); end
  else
      dx = 1;
  end

  % Y spacing
  y = Y(:,1); y = y(isfinite(y));
  if numel(y) >= 2
      dyv = diff(y); dyv = dyv(isfinite(dyv) & dyv~=0);
      if isempty(dyv), dy = 1; else, dy = median(dyv); end
  else
      dy = 1;
  end

  if ~isfinite(dx) || dx<=0, dx = 1; end
  if ~isfinite(dy) || dy<=0, dy = 1; end

  % Use min/max so extents are correct whether grids ascend or descend
  xExt = [min(x)-dx/2, max(x)+dx/2];
  yExt = [min(y)-dy/2, max(y)+dy/2];
end

% Restricted Cubic Spline basis
function X = rcs_basis(x, knots)
  x = x(:); k = knots(:).'; K = numel(k);
  if K < 4, error('Need at least 4 knots for RCS'); end
  X = zeros(numel(x), K-2);
  pp = @(z) max(z,0).^3;
  denom = (k(K) - k(K-1));
  for j = 1:(K-2)
    X(:,j) = pp(x - k(j)) ...
           - pp(x - k(K-1)) * ((k(K) - k(j))/denom) ...
           + pp(x - k(K))   * ((k(K-1) - k(j))/denom);
  end
end

% Robust 2–98% limits ignoring NaNs (optional mask)
function [lo, hi] = robust_clim_fast(Z, M)
  if nargin < 2 || isempty(M)
    Zs = Z(isfinite(Z));
  else
    Zs = Z(isfinite(Z) & M);
  end
  if isempty(Zs), lo = 0; hi = 1; return; end
  Zs = double(Zs);
  lo = prctile(Zs, 2); hi = prctile(Zs, 98);
end

% Simple diverging colormap (blue–white–red)
function C = diverging_cmap(n)
  if nargin<1, n=256; end
  m = floor(n/2);
  up = [linspace(0,1,m)' linspace(0,1,m)' ones(m,1)];   % blue -> white
  dn = [ones(n-m,1) linspace(1,0,n-m)' linspace(1,0,n-m)']; % white -> red
  C  = [up; dn];
end

% Pretty 2D curve with per-vertex color
function plot_gradient_curve(ax, x, y, c, lw)
  x = x(:)'; y = y(:)'; c = c(:)';
  if numel(x) < 2, return; end
  X = [x; x]; Y = [y; y]; Z = zeros(2, numel(x)); C = [c; c];
  surface(ax, X, Y, Z, C, 'FaceColor','none', 'EdgeColor','interp', 'LineWidth', lw);
end

% MD5 of a logical array (size + data)
function h = md5_of_logical(L)
    if isempty(L), h = repmat('0',1,32); return; end
    bs   = getByteStreamFromArray(logical(L));
    md   = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(bs));
    raw  = typecast(md.digest(),'uint8');
    h    = lower(reshape(dec2hex(raw,2).',1,[]));
end

% Pr[(M - Gmin) > 0] when (M - Gmin) ~ N(mu, sig_m^2 + sig_g^2)
function pwet = wet_prob_from_sigmas(gdif_vec, sig_m, sig_g)
  if isempty(sig_g), sig_g = 0; end
  s2  = sig_m.^2 + sig_g.^2;
  s   = sqrt(max(s2, eps));
  pwet = 0.5 * (1 + erf((gdif_vec ./ s) ./ sqrt(2)));
end

% Finite-difference uncertainty propagation for Gmin
% sigma(Gmin)^2 ≈ (∂G/∂Ts)^2 σ_Ts^2 + (∂G/∂Mb)^2 σ_Mb^2 + (∂G/∂H)^2 σ_H^2
function [sigma_map, sigma_scalar] = estimate_sigma_gmin_fd(S, cfg)
  
    key = sprintf('%s_%s_%s', md5_of_logical(isfinite(S.H)), ...
                         md5_of_logical(S.Ts>-Inf), md5_of_logical(S.Mb>-Inf));
    cache_path = fullfile(cfg.cache_dir, ['sigma_gmin_fd_' key '.mat']);
    if isfile(cache_path)
        L = load(cache_path, 'sigma_map', 'sigma_scalar');
        sigma_map = L.sigma_map; sigma_scalar = L.sigma_scalar; return;
    end

  % Baseline Gmin
  [~, G0, ~] = processGHF(S.models.Hazzard, S.Ts, S.Mb, S.H); %#ok<ASGLU>

  % Ts derivative
  fprintf('   Ts +/- %g ...\n', cfg.gmin_unc.delta_Ts);
  dT = cfg.gmin_unc.delta_Ts;
  [~, Gp, ~] = processGHF(S.models.Hazzard, S.Ts + dT, S.Mb, S.H);
  [~, Gm, ~] = processGHF(S.models.Hazzard, S.Ts - dT, S.Mb, S.H);
  dGdT = (Gp - Gm) / (2*dT);

  % Mb derivative
  fprintf('   Mb +/- %g ...\n', cfg.gmin_unc.delta_Mb);
  dA = cfg.gmin_unc.delta_Mb;
  [~, Gp, ~] = processGHF(S.models.Hazzard, S.Ts, S.Mb + dA, S.H);
  [~, Gm, ~] = processGHF(S.models.Hazzard, S.Ts, S.Mb - dA, S.H);
  dGdA = (Gp - Gm) / (2*dA);

  % H derivative
  fprintf('   H +/- %g ...\n', cfg.gmin_unc.delta_H);
  dH = cfg.gmin_unc.delta_H;
  [~, Gp, ~] = processGHF(S.models.Hazzard, S.Ts, S.Mb, S.H + dH);
  [~, Gm, ~] = processGHF(S.models.Hazzard, S.Ts, S.Mb, S.H - dH);
  dGdH = (Gp - Gm) / (2*dH);

  % Combine in quadrature
  sT = cfg.gmin_unc.sigma_Ts;
  sA = cfg.gmin_unc.sigma_Mb;
  sH = cfg.gmin_unc.sigma_H;

  sigma_map = sqrt( (dGdT .* sT).^2 + (dGdA .* sA).^2 + (dGdH .* sH).^2 );

  % ROI RMS scalar summary
  if isfield(S,'roi_mask') && ~isempty(S.roi_mask)
      vv = sigma_map(S.roi_mask & isfinite(sigma_map));
  else
      vv = sigma_map(isfinite(sigma_map));
  end
  sigma_scalar = sqrt(mean(double(vv).^2, 'omitnan'));

  try save(cache_path, 'sigma_map','sigma_scalar','-v7.3'); catch
  end
end

function val = get_field_default(s, fld, defaultVal)
    if isstruct(s) && isfield(s, fld) && ~isempty(s.(fld))
        val = s.(fld);
    else
        val = defaultVal;
    end
end

function assert_in_range(x, lo, hi, name)
    bad = ~isnan(x) & (x<lo | x>hi);
    assert(~any(bad(:)), '%s out of range [%g,%g] (violations=%d)', name, lo, hi, nnz(bad));
end

function assert_all_finite(x, name)
    assert(all(isfinite(x(:)) | isnan(x(:))), '%s contains +/-Inf (not allowed)', name);
end

function math_warn(cfg, msg, varargin)
    dbg(1,cfg,['[MATH WARN] ' msg], varargin{:});
end

function ok = check_confusion_identities(m)
    N2 = m.TP + m.FP + m.TN + m.FN;
    ok = true;
    ok = ok & (abs(m.N - N2) <= 1e-9);
    v = [m.ACC m.PREC m.SPEC m.REC m.F1];
    ok = ok & all(v >= -1e-9 & v <= 1+1e-9);
end

function fig = plot_overlay_on_H(X_km, Y_km, H, Z, climZ, mode, cfg, mask_slow, varargin)
% plot_overlay_on_H - H as grayscale background, Z overlaid with alpha
% Inputs are unchanged from your version.

    p = inputParser;
    addParameter(p,'SinkMask',[],@(x)islogical(x));
    addParameter(p,'SinkStyle',struct(),@isstruct);
    parse(p,varargin{:});
    sinkMask  = p.Results.SinkMask;

    % --- extents so all rasters share the same placement ---
    [xExt, yExt] = grid_extents(X_km, Y_km);

    % --- pick colormap for Z (H will use gray) ---
    if strcmpi(mode,'diverging')
        cmapZ = diverging_cmap(256);
    else
        cmapZ = parula(256);
    end

    % --- figure + two overlaid axes (separate colormaps) ---
    fig = figure('Visible', 'off', 'Color','w','Position',[100 100 900 800]);

    % Bottom axis: H grayscale background
    axH = axes('Parent',fig); hold(axH,'on'); axis(axH,'image'); box(axH,'on');
    set(axH,'YDir','normal','Layer','bottom');
    % robust 2–98% stretch for H
    Hfinite = H(isfinite(H));
    if isempty(Hfinite)
        Hlo = 0; Hhi = 1;
    else
        Hlo = prctile(double(Hfinite),2);
        Hhi = prctile(double(Hfinite),98);
    end
    Himg = mat2gray(double(H), [Hlo Hhi]);
    imagesc(axH, xExt, yExt, Himg);
    colormap(axH, gray);  % grayscale for H
    % no colorbar on H axis

    % Top axis: Z with alpha over the H background
    axZ = axes('Parent',fig, 'Color','none'); hold(axZ,'on'); axis(axZ,'image'); box(axZ,'on');
    set(axZ,'YDir','normal','Layer','top');
    him = imagesc(axZ, xExt, yExt, Z);
    % alpha from mask_slow
    if ~isempty(mask_slow)
        set(him,'AlphaData', double(mask_slow) * cfg.overlay_alpha, 'AlphaDataMapping','none');
    else
        set(him,'AlphaData', cfg.overlay_alpha);
    end
    colormap(axZ, cmapZ);
    clim(axZ, climZ);
    cb = colorbar(axZ); ylabel(cb,'mW m^{-2}');

    % --- keep only ONE set of ticks/labels/box (top axis) ---
    set(axH, 'XTick',[], 'YTick',[], 'XColor','none', 'YColor','none', 'Box','off');  % bottom axis silent
    set(axZ, 'Box','on');  % show the box only once (top axis)
    
    % --- hard-align the axes rectangles ---
    set(axZ, 'Units','normalized');
    set(axH, 'Units','normalized');
    set(axZ, 'Position', get(axH,'Position'));
    
    % --- and lock limits together to avoid drift ---
    if isgraphics(axH,'axes') && isgraphics(axZ,'axes')
        lh = linkprop([axH, axZ], {'XLim','YLim'});
        % Keep the link object alive for the life of the figure:
        setappdata(fig, 'linkXY', lh);
    else
        warning('plot_overlay_on_H:axesInvalid','Skipping axis linking: an axis handle is invalid.');
    end
    axis(axZ,'image');   % keep equal aspect on the top axis

    % H contours on the top axis (so lines sit above the Z image)
    if isfield(cfg,'H_contour_interval') && ~isempty(cfg.H_contour_interval)
        try
            [C,hC] = contour(axZ, X_km, Y_km, H, ...
                             'LineColor',[0.3 0.3 0.3], 'LineWidth',0.5);
            clabel(C,hC,'LabelSpacing',400,'FontSize',8,'Color',[0.3 0.3 0.3]);
        catch
        end
    end

    % ---- ROI / sink mask as green fill + outlines on TOP axis ----
    if ~isempty(sinkMask) && any(sinkMask(:))
        % draw using the SAME X/Y vectors as the background, not [xExt yExt]
        RGB = cat(3, 0.10*ones(size(sinkMask)), ...
                     0.65*ones(size(sinkMask)), ...
                     0.20*ones(size(sinkMask)));
        hmask = imagesc(axZ, X_km(1,:), Y_km(:,1), RGB);  % <— imagesc + grid vectors
        set(hmask, 'AlphaData', 0.40 * double(sinkMask), ...   % a bit stronger
                   'AlphaDataMapping','none', ...
                   'HitTest','off', 'PickableParts','none');
    
        uistack(hmask, 'top'); 
    
        % crisp outlines 
        B = bwboundaries(sinkMask, 8);
        for b = 1:numel(B)
            bb = B{b};
            plot(axZ, X_km(1,bb(:,2)), Y_km(bb(:,1),1), '-', ...
                 'Color',[0.10 0.65 0.20], 'LineWidth',1.2, ...
                 'HitTest','off');
        end
    else
        fprintf('SinkMask empty or all-false; nothing to draw.\n');
    end

    xlabel(axZ,'X (km)'); ylabel(axZ,'Y (km)');
    title(axZ,'Overlay on H (grayscale background)');
end

function m = metrics_from_counts(TP,FP,TN,FN)
% Robust metrics from (expected) counts; safe for tiny denominators
    N = TP + TN + FP + FN + eps;
    PREC = TP / max(TP + FP, eps);
    SPEC = TN / max(TN + FP, eps);
    REC  = TP / max(TP + FN, eps);
    ACC  = (TP + TN) / N;
    F1   = 2 * PREC * REC / max(PREC + REC, eps);
    m = struct('PREC',PREC,'SPEC',SPEC,'REC',REC,'ACC',ACC,'F1',F1);
end

function lim = stat_axis_limits(stat_name)
% Axis limits per metric
    s = upper(stat_name);
    switch s
        case 'MCC'
            lim = [-1 1];
        otherwise
            lim = [0 1];
    end
end

function plot_colored_ci_pair(S, cfg, statX, statY, mode)
% plot colored CI rectangles for (statX, statY)
% mode = 'adj' or 'raw' ; stat names like 'SPEC','REC','PREC','ACC','F1'

    statX = upper(statX); statY = upper(statY); mode = lower(mode);
    GX = S.Gstats.([statX '_' mode]);
    GY = S.Gstats.([statY '_' mode]);

    X  = GX.EXP;   Xlohi  = GX.CI;   Xlohi_i = GX.CIi;
    Y  = GY.EXP;   Ylohi  = GY.CI;   Ylohi_i = GY.CIi;

    nM = numel(S.names);
    is_flat = startsWith(S.names,'Flat_');
    colors = lines(max(1,nM));
    epsw = 1e-6;

    % Figure/axes
    f = figure('Color','w','Position',[120 120 900 820]);
    ax = axes(f); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); ax.Layer='top';
    xlabel(ax, sprintf('%s (%s, expected)', statX, upper(mode)));
    ylabel(ax, sprintf('%s (%s, expected)', statY, upper(mode)));
    xlim(ax, stat_axis_limits(statX));
    ylim(ax, stat_axis_limits(statY));

    % ---- Flat GHF sweep (gradient by flat value in mW m^-2) ----
    If = find(is_flat);
    xs = X(If); ys = Y(If); vals = S.model_cvals(If);  % flat GHF values
    good = isfinite(xs) & isfinite(ys) & isfinite(vals);
    xs = xs(good); ys = ys(good); vals = vals(good);

    if numel(vals) >= 2
        % sort by flat value so the gradient is monotonic
        [vals_s, ord] = sort(vals);
        xs_s = xs(ord); ys_s = ys(ord);

        if isfield(cfg.flat,'interp_enable') && cfg.flat.interp_enable && numel(vals_s) >= 2
            % interpolate the curve parameterized by "flat value"
            npts = max(20, cfg.flat.interp_points);
            [u, iu] = unique(vals_s);                  % ensure strictly increasing param
            vi = linspace(min(u), max(u), npts);       % parameter = flat value (mW m^-2)
            mth = 'spline';
            if isfield(cfg.flat,'interp_method') && ~isempty(cfg.flat.interp_method)
                mth = cfg.flat.interp_method;
            end
            xi = interp1(u, xs_s(iu), vi, mth, 'extrap');
            yi = interp1(u, ys_s(iu), vi, mth, 'extrap');

            % gradient line using your helper (colored by flat value)
            plot_gradient_curve(ax, xi, yi, vi, 15);
        else
            % no interpolation; color per-vertex by its flat value
            plot_gradient_curve(ax, xs_s, ys_s, vals_s, 15);
        end

        % set colormap & color scale for the gradient (others use explicit RGB)
        colormap(ax, diverging_cmap(256));
        clim(ax, [min(vals) max(vals)]);
        grid off;
        cb = colorbar(ax);
        ylabel(cb, 'Flat GHF (mW m^{-2})');
    end

        % Non-flat models: colored CI envelopes + markers + labels
    for i = 1:nM
        if is_flat(i), continue; end
        if ~all(isfinite([X(i) Y(i)])), continue; end

        % outer band
        if all(isfinite([Xlohi(i,:) Ylohi(i,:)]))
            x1 = Xlohi(i,1); x2 = Xlohi(i,2); if x2<=x1, x2 = x1+epsw; end
            y1 = Ylohi(i,1); y2 = Ylohi(i,2); if y2<=y1, y2 = y1+epsw; end
            patch(ax, [x1 x2 x2 x1], [y1 y1 y2 y2], colors(i,:), ...
                  'FaceAlpha',0.12, 'EdgeColor','none', 'HitTest','off');
        end

        % inner band
        if all(isfinite([Xlohi_i(i,:) Ylohi_i(i,:)]))
            x1 = Xlohi_i(i,1); x2 = Xlohi_i(i,2); if x2<=x1, x2 = x1+epsw; end
            y1 = Ylohi_i(i,1); y2 = Ylohi_i(i,2); if y2<=y1, y2 = y1+epsw; end
            patch(ax, [x1 x2 x2 x1], [y1 y1 y2 y2], colors(i,:), ...
                  'FaceAlpha',0.25, 'EdgeColor','none', 'HitTest','off');
        end

        % point + label
        plot(ax, X(i), Y(i), 'o', 'MarkerFaceColor',colors(i,:), ...
             'MarkerEdgeColor','k', 'MarkerSize',6, 'HandleVisibility','off');
        if i <= numel(S.titles)
            dx = 0.01 * range(xlim(ax)); dy = 0.01 * range(ylim(ax));
            text(ax, X(i)+dx, Y(i)+dy, S.titles{i}, ...
                 'FontSize',10, 'Color',colors(i,:), 'Interpreter','none');
        end
    end

    title(ax, sprintf('%s – %s (%s): colored CIs (models) + flat sweep', ...
        statX, statY, upper(mode)));
end

function t = dbg_tic(label, cfg)
    dbg(1,cfg,'[TIMER START] %s', label); t = tic;
end
function dbg_toc(t, label, cfg)
    dbg(1,cfg,'[TIMER END]   %s — %.3f s', label, toc(t));
end

function dbg_progress(i, n, cfg, every)
    if nargin<4 || isempty(every), every = 25; end
    if i==1 || i==n || mod(i,every)==0
        dbg(2,cfg,'Progress: %d/%d (%.1f%%)', i, n, 100*i/n);
        drawnow;
    end
end

function dbg_pause(cfg, why)
    persistent skip_all; if isempty(skip_all), skip_all = false; end
    if ~isfield(cfg,'debug') || ~cfg.debug.enable, return; end
    if isfield(cfg.debug,'pause_on_fig') && ~cfg.debug.pause_on_fig, return; end
    if skip_all, return; end
    if nargin<2, why = 'debug'; end
    resp = input(sprintf('Paused after %s. [Enter]=continue, s=skip all, q=abort > ', why),'s');
    if strcmpi(resp,'q'), error('User aborted run.'); end
    if strcmpi(resp,'s'), skip_all = true; end
end

function dbg_checkpoint(label, varsStruct, cfg)
    if ~isfield(cfg,'debug') || ~cfg.debug.enable, return; end
    if ~isfield(cfg.debug,'save_snapshots') || ~cfg.debug.save_snapshots, return; end
    if ~isfield(cfg.debug,'snapshot_dir') || isempty(cfg.debug.snapshot_dir)
        cfg.debug.snapshot_dir = fullfile(cfg.outdir,'debug');
    end
    if ~exist(cfg.debug.snapshot_dir,'dir'), mkdir(cfg.debug.snapshot_dir); end
    fname = fullfile(cfg.debug.snapshot_dir, sprintf('snap_%s_%s.mat', ...
            datetime('now','yyyymmdd_HHMMSS'), regexprep(label,'\W+','_')));
    try
        save(fname, '-struct', 'varsStruct','-v7');
        dbg(2,cfg,'[CHECKPOINT] %s -> %s', label, fname);
    catch ME
        dbg(1,cfg,'[CHECKPOINT FAILED] %s: %s', label, ME.message);
    end
end

function quicklook_field(name, X_km, Y_km, F, mask, cfg)
% quicklook_field - diagnostic plot for a raster F on the grid
%   name      : title string
%   X_km,Y_km : grid (km)
%   F         : field to show (same size as grid)
%   mask      : optional logical mask to overlay (e.g., sinks/ROI)
%   cfg       : config (uses cfg.debug.quicklook_pause)

    if ~isfield(cfg,'debug') || ~cfg.debug.enable || ~isfield(cfg.debug,'quicklook_all') || ~cfg.debug.quicklook_all
        return;
    end
    if isempty(F) || ~isequal(size(F), size(X_km)) || ~isequal(size(F), size(Y_km))
        return;
    end

    ax = nexttile;                 

    imagesc(ax, X_km(1,:), Y_km(:,1), F);
    set(ax,'YDir','normal'); axis(ax,'image');
    title(ax, name, 'Interpreter','none');
    colorbar(ax); colormap(ax, parula);

    if nargin>=5 && ~isempty(mask)
        [sy,sx] = find(mask);
        hold(ax,'on');
        % style from cfg.sinks if available
        sty = struct('marker','.', 'size',4, 'color',[0.05 0.05 0.05], 'alpha',0.85);
        if isfield(cfg,'sinks')
            if isfield(cfg.sinks,'marker'), sty.marker = cfg.sinks.marker; end
            if isfield(cfg.sinks,'size'),   sty.size   = cfg.sinks.size;   end
            if isfield(cfg.sinks,'color'),  sty.color  = cfg.sinks.color;  end
            if isfield(cfg.sinks,'alpha'),  sty.alpha  = cfg.sinks.alpha;  end
        end
        scatter(ax, X_km(1,sx), Y_km(sy,1), sty.size, sty.marker, ...
                'MarkerEdgeColor',sty.color, ...
                'MarkerEdgeAlpha',sty.alpha, 'MarkerFaceAlpha',sty.alpha);
    end

    if isfield(cfg.debug,'quicklook_pause') && cfg.debug.quicklook_pause
        dbg_pause(cfg, ['quicklook: ' name]);
    end
end

function o = override_struct(base, over)
% override_struct - shallow merge: fields in 'over' override 'base'
    o = base;
    if nargin<2 || isempty(over) || ~isstruct(over), return; end
    f = fieldnames(over);
    for ii = 1:numel(f)
        o.(f{ii}) = over.(f{ii});
    end
end