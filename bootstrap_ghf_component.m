%% ================================================================
%  GHF model evaluation
%
%  Inputs expected:
%    datasets_for_gmin/GHF_*_interp.mat              (GHF fields)
%    datasets_for_gmin/UNC_*_interp.mat              (σ fields, where applicable)
%    datasets_for_gmin/BMIN_Losing_interp.mat        (min GHF, Losing)
%    datasets_for_gmin/BMAX_Losing_interp.mat        (max GHF, Losing)
%    datasets_for_gmin/specularity.mat               (Q)
%    datasets_for_gmin/coldex_icethk.mat             (Xgrid, Ygrid, H)
%    datasets_for_gmin/Ts_interp.mat                 (Ts_interp)
%    datasets_for_gmin/Mb_interp.mat                 (Mb_interp)
%    datasets_for_gmin/mouginot_icevel.mat           (icevel)
%    datasets_for_gmin/sink_mask_new.mat (logical)   var: sink_mask
%  ================================================================
function S = bootstrap_ghf_component()

check_dependencies(struct('autoInstall',true,'needAMT',false,'verbose',true));
addpath('datasets_for_gmin');
safe = @(s) regexprep(char(s), '[^A-Za-z0-9._-]+', '_');  % sanitize for filenames

%% -------------------- CONFIG --------------------
cfg = struct( ...
  'spec_thresh',        0.20, ...   % Specularity content threshold
  'v_keep',             3.5, ...    % Ice velocity threshold (m/yr)
  'target_pos_pct',     15, ...     % Percent of observations to accept as "wet" after icethk bias
  'decision_margin',    0.1, ...    % Value of Gdif-Gmin above or below 0 to accept as "wet"
  'outdir',             'figs_out', ...
  'overwrite',          true, ...
  'overlay_alpha',      0.4, ...
  'font_size',          20, ...
  'to_single',          true ...
);

% === Region split from boundary polyline (ISPB/OSPB/all) ===
cfg.region = struct( ...
  'mode',          'all', ...         % 'ISPB' | 'OSPB' | 'all'
  'left_label',    'ISPB', ...         % grid-left (toward x=0)git
  'right_label',   'OSPB', ...
  'cache',         true ...
);

% Flat GHF sweep
cfg.flat = struct( ...
  'enable',         true, ...
  'values',         20:1:80, ...        % mW/m^2
  'interp_enable',  true, ...
  'interp_points',  200, ...
  'interp_method',  'spline' ...
);

% Sinks (used to define ROI; WITH connectivity grouping)
cfg.sinks = struct( ...
  'marker',        '.', ...
  'size',          4, ...
  'color',         [0.05 0.05 0.05], ...
  'alpha',         0.85, ...
  'min_size',      10, ...              % Number of pixels
  'connectivity',  8 ...
);

% Uncertainty — analytic expectation + sink bootstrap
cfg.uncertainty = struct( ...
  'mode',        'analytic', ...        % (reserved)
  'bootstrap_mode', 'component', ...    % 'pixel' | 'component'
  'n_boot',      500, ...               % Number of component resamples
  'seed',        42, ...
  'mfrac',       0.8, ...               % Fraction of components to resample
  'band_main',   0.95, ...              % outer shaded band (e.g., 90%)
  'band_inner',  0.50 ...               % inner shaded band (e.g., 50%)
);

% Gmin uncertainty via finite-difference propagation (per pixel)
cfg.gmin_unc = struct( ...
  'enable',     true, ... 
  'sigma_Ts',   1.5, ...                % Surface temperature (K) van Dalum 
  'sigma_Mb',   0.02, ...               % Mass balance (m/yr IE) van Dalum 
  'sigma_H',    50.0, ...               % Ice thickness (m)
  'delta_Ts',   0.5, ...
  'delta_Mb',   0.01, ...
  'delta_H',    5.0 ...
);

fallback_err = struct('FoxMaule',21,'Hazzard',8,'Martos',10,'Shen',9,'Stal',12,'An',45,'Losing',25);  % σ (mW/m^2)

% Global rectangle mask in km
cfg.rect = struct( ...
  'enable',  true, ...                 % true to activate
  'xlim',    [10 700], ...             % [xmin xmax] alternative to 'rects'
  'ylim',    [-200 350], ...           % [ymin ymax]
  'rects',   [], ...                   % Nx4 rows: [xmin xmax ymin ymax];
  'mode',    'include' ...             % 'include' | 'exclude'
);

cfg.plot_labels = struct( ...
  'avoid_overlap', true, ...   % enable repel
  'max_iter',      200, ...    % iterations for nudging
  'step',          0.01, ...   % fraction of axis range per nudge
  'pad',           0.005 ...   % extra padding between labels (axis frac)
);

if ~exist(cfg.outdir,'dir'), mkdir(cfg.outdir); end

%% -------------------- LOAD -------------------- %%
fprintf('Loading data ...\n');
S = struct();

% GHF models
t=load('datasets_for_gmin/GHF_Hazzard_interp.mat','GHF_Hazzard_interp');  GHF_Hazzard_interp=t.GHF_Hazzard_interp; clear t
t=load('datasets_for_gmin/GHF_Martos_interp.mat','GHF_Martos_interp');    GHF_Martos_interp=t.GHF_Martos_interp;   clear t
t=load('datasets_for_gmin/GHF_Shen_interp.mat','GHF_Shen_interp');        GHF_Shen_interp=t.GHF_Shen_interp;       clear t
t=load('datasets_for_gmin/GHF_Stal_interp.mat','GHF_Stal_interp');        GHF_Stal_interp=t.GHF_Stal_interp;       clear t
t=load('datasets_for_gmin/GHF_Losing_interp.mat','GHF_Losing_interp');    GHF_Losing_interp=t.GHF_Losing_interp;   clear t
t=load('datasets_for_gmin/GHF_An_interp.mat','GHF_An_interp');            GHF_An_interp=t.GHF_An_interp;           clear t
t=load('datasets_for_gmin/GHF_FoxMaule_interp.mat','GHF_FoxMaule_interp');GHF_FoxMaule_interp=t.GHF_FoxMaule_interp;clear t

S.models = struct('FoxMaule',GHF_FoxMaule_interp,'Hazzard',GHF_Hazzard_interp, ...
    'Martos',GHF_Martos_interp,'Shen',GHF_Shen_interp,'Stal', ...
    GHF_Stal_interp,'An',GHF_An_interp,'Losing',GHF_Losing_interp);
S.names  = {'FoxMaule','Hazzard','Martos','Shen','Stal','An','Losing'};
S.titles = {'Fox Maule','Hazzard','Martos','Shen','Stål','An','Lösing'};

t = load('datasets_for_gmin/coldex_icethk.mat','Xgrid','Ygrid','H');  % Ice thickness + grid
S.Xgrid = t.Xgrid; S.Ygrid = t.Ygrid; S.H = t.H; clear t 
S.X_km = S.Xgrid/1000; S.Y_km = S.Ygrid/1000;
% S.H = S.h - 200

t = load('datasets_for_gmin/specularity.mat','Q'); S.Q = t.Q; clear t  % Specularity

t = load('datasets_for_gmin/Ts_interp.mat','Ts_interp'); S.Ts = t.Ts_interp; clear t  % Surface temp

t=load('datasets_for_gmin/Mb_interp.mat','Mb_interp'); S.Mb=t.Mb_interp; clear t  % Mass balance
% S.Mb = S.Mb * 0.6;

t = load('datasets_for_gmin/sink_mask_new.mat', 'sink_mask'); S.sink_mask = logical(t.sink_mask); clear t  % Sinks

t=load('datasets_for_gmin/mouginot_icevel.mat','icevel'); S.icevel=t.icevel; clear t  % Ice Velocity

% === Build ISPB/OSPB split from polyline ===
B = shaperead('datasets_for_gmin/roughness_boundary.shp');   
bx = [B.X]'; by = [B.Y]';
good = isfinite(bx) & isfinite(by);
bx = bx(good) / 1000;    % km
by = by(good) / 1000;    % km

[S.mask_ISPB, S.mask_OSPB, ~] = build_spb_split_masks(S.X_km, S.Y_km, [bx by]);
fprintf('[ISPB] pixels=%d  [OSPB] pixels=%d\n', nnz(S.mask_ISPB), nnz(S.mask_OSPB));

S.rect_mask = build_rect_mask(S.X_km, S.Y_km, cfg.rect);

xy_mask = (S.mask_ISPB | S.mask_OSPB) & S.rect_mask;

% Choose region
switch lower(cfg.region.mode)
  case 'ispb', REG_MASK = S.mask_ISPB & S.rect_mask;
  case 'ospb', REG_MASK = S.mask_OSPB & S.rect_mask;
  case 'all',  REG_MASK = xy_mask;
  otherwise, error('cfg.region.mode must be ISPB | OSPB | all');
end

% --- Per-pixel mean of all non-flat (original) models ---
nonflat_names = {'FoxMaule','Hazzard','Martos','Shen','Stal','An','Losing'};

% stack and mean (omit NaNs)
stk = cat(3, S.models.(nonflat_names{1}), S.models.(nonflat_names{2}), ...
         S.models.(nonflat_names{3}), S.models.(nonflat_names{4}), ...
         S.models.(nonflat_names{5}), S.models.(nonflat_names{6}), ...
         S.models.(nonflat_names{7}));
S.models.MeanNonFlat = mean(stk, 3, 'omitnan');

% register like every other model
S.names{end+1}  = 'MeanNonFlat';
S.titles{end+1} = 'Mean (non-flat)';

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

% --- Uncertainty for MeanNonFlat ---
unc_fields = {'Hazzard','Martos','Shen','Stal'};  % models that ship with UNC_* grids
have_unc = unc_fields(isfield(S.unc, unc_fields));
if ~isempty(have_unc)
  Ustk = [];
  for k = 1:numel(have_unc)
      if ~isempty(S.unc.(have_unc{k})), Ustk = cat(3, Ustk, S.unc.(have_unc{k})); end
  end
  if ~isempty(Ustk), S.unc.MeanNonFlat = mean(Ustk, 3, 'omitnan'); end
end

% If no grid was created above, provide a scalar fallback for the mean
if ~isfield(S.unc,'MeanNonFlat') || isempty(S.unc.MeanNonFlat)
  mfallback = mean([fallback_err.FoxMaule, fallback_err.Hazzard, fallback_err.Martos, ...
                  fallback_err.Shen, fallback_err.Stal,   fallback_err.An, ...
                  fallback_err.Losing], 'omitnan');
  fallback_err.MeanNonFlat = mfallback;   % picked up by get_field_default(...) later
end

% Cast → single 
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

%% -------------------- ROI MASK (ISPB/OSPB/all + slow flow + sinks) -----
% Validity for upstream fields (you can add more if needed)
valid_mask = isfinite(S.Q) & isfinite(S.H) & isfinite(S.icevel);

% Slow flow (keep <= v_keep)
fast_mask  = isfinite(S.icevel) & (S.icevel > cfg.v_keep);
slow_mask  = valid_mask & ~fast_mask;

% Final ROI: valid ∧ slow ∧ chosen region ∧ sinks (sinks are mandatory)
S.roi_mask = slow_mask & REG_MASK & S.sink_mask;

% Visualization mask for maps (NO sinks): valid ∧ slow ∧ chosen region ∧ rect
S.viz_mask = slow_mask & REG_MASK;

% Keep ROI and EVAL masks identical
idx = find(S.roi_mask);
S.eval_mask_full = S.roi_mask;    % grid-space logical
S.eval_mask      = true(numel(idx),1);   % ROI-vector space is simply "all true"

% Labels use exactly ROI, same as everything downstream
H_vec     = S.H(idx);
y_raw_vec = (S.Q(idx) > cfg.spec_thresh);   % unchanged

% --- Connected components used for bootstrap (true sink components only) ---
comp_mask = S.roi_mask & S.sink_mask;  % only true sink pixels inside ROI
conn = cfg.sinks.connectivity;         % 4 or 8
CC = bwconncomp(comp_mask, cfg.sinks.connectivity);

% prune tiny components (only when using sinks)
keep = cellfun(@numel, CC.PixelIdxList) >= cfg.sinks.min_size;
if ~all(keep)
    comp_mask(vertcat(CC.PixelIdxList{~keep})) = false;
    CC = bwconncomp(comp_mask, conn);
end

S.comp_id = zeros(size(S.roi_mask),'uint32');
for k = 1:CC.NumObjects, S.comp_id(CC.PixelIdxList{k}) = k; end
S.n_sinks = CC.NumObjects;

% Cache dir keyed to ROI + sinks fingerprint (needed below)
roi_fp  = md5_of_logical(S.roi_mask);
sink_fp = md5_of_logical(S.sink_mask);
cfg.cache_dir = fullfile(cfg.outdir, sprintf('roi_cache_%s_%s', roi_fp(1:12), sink_fp(1:12)));
if ~exist(cfg.cache_dir,'dir'), mkdir(cfg.cache_dir); end

%% -------------------- LABELS (RAW + THK-ADJ) --------------------
fprintf('Fitting RCS GLM and computing z-residuals...\n');
Hz   = zscore(log(double(H_vec) + 1));
q    = prctile(Hz, [5 27.5 50 72.5 95]);
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

%% -------------------- PROCESS GHF → ROI cache --------------------
nM = numel(S.names);
gdif_cache = cell(nM,1);
S.model_cvals = nan(nM,1);

for i = 1:nM
fprintf('[Uncertainty] Model %d/%d: %s\n', i, nM, S.names{i});
name_i = S.names{i};
Mi = S.models.(name_i);
if startsWith(name_i,'Flat_'), v = sscanf(name_i,'Flat_%f'); S.model_cvals(i) = v;
    else, S.model_cvals(i) = median(Mi(idx),'omitnan');
end
end

fprintf('Processing models to Gdif (ROI cache)...\n');

% Gmin (one-time)
gmin_path = fullfile(cfg.cache_dir, 'Gmin_full.mat');
[~, Gmin_full, ~] = processGHF(S.models.Hazzard, S.Ts, S.Mb, S.H); % model choice irrelevant for Gmin
    try save(gmin_path, 'Gmin_full', '-v7.3'); catch, warning('Could not cache Gmin_full'); end
gmin_vec = Gmin_full(idx);

% Per-model GΔ (ROI vectors) + cache
for i = 1:nM
    cache_i = fullfile(cfg.cache_dir, sprintf('gdifvec_%s.mat', S.names{i}));
    if isfile(cache_i)
        L = load(cache_i, 'gdif_vec'); gdif_cache{i} = L.gdif_vec;
    else
        Mi = S.models.(S.names{i});
        try
            [~, ~, Gdif_full] = processGHF(Mi, S.Ts, S.Mb, S.H);
        catch ME 
            disp(getReport(ME,'extended')); rethrow(ME);
        end
        gdif_vec = Gdif_full(idx);
        gdif_cache{i} = gdif_vec;
        try save(cache_i, 'gdif_vec', '-v7.3'); catch, end
    end
end

% Keep ROI and EVAL identical by design
S.eval_mask       = true(numel(idx),1);
S.eval_mask_full  = S.roi_mask;

% ---- MATH CHECKS: eval mask coherence ----
M_eval = S.eval_mask;
assert(nnz(M_eval)>0, 'eval_mask is empty');
assert(numel(M_eval)==numel(gmin_vec), 'eval_mask length mismatch');

for i=1:nM, assert(numel(gdif_cache{i})==numel(gmin_vec), 'gdif_cache len mismatch'); end

%% === Bootstrap groups (connected sinks) precompute (for 'component' mode) ===
S.boot_groups = struct('nC',0,'groups',{{}},'uniqC',[]);
if strcmpi(cfg.uncertainty.bootstrap_mode,'component')
    Mglob = S.eval_mask;
    compM = S.comp_id(idx);           % component id per ROI vector pixel
    compM = compM(Mglob);             % restrict to eval-mask space
    if any(compM)
        [uniqC,~,gidx] = unique(compM(compM>0));
        groups = accumarray(gidx, find(compM>0), [], @(v){v});
        S.boot_groups.uniqC  = uniqC;
        S.boot_groups.groups = groups;
        S.boot_groups.nC     = numel(groups);
    end
end

%% -------------------- Estimate σ(Gmin) via finite differences -------------
if cfg.gmin_unc.enable 
    [S.sigma_gmin, S.sigma_gmin_scalar] = estimate_sigma_gmin_fd(S, cfg);
else
    S.sigma_gmin = [];
    S.sigma_gmin_scalar = NaN;
end

% ---- MATH CHECKS: sigma_gmin ----
if ~isempty(S.sigma_gmin), assert(all(S.sigma_gmin(:) >= 0 | isnan(S.sigma_gmin(:))), 'sigma_gmin must be >= 0'); end

%% -------------------- EVALUATE (point metrics, no uncertainty) ------------
fprintf('Evaluating point metrics on ROI...\n');

TP_r  = zeros(nM,1,'double'); FP_r  = TP_r; TN_r = TP_r; FN_r = TP_r;
TP_a  = TP_r; FP_a  = TP_r;   TN_a  = TP_r; FN_a = TP_r;
ACC_r = TP_r; ACC_a = TP_r;   PR_r  = TP_r; PR_a = TP_r; RC_r = TP_r; RC_a = TP_r; F1_r = TP_r; F1_a = TP_r;
TPFP_r= TP_r; TNFN_r= TP_r;   TPFP_a= TP_r; TNFN_a= TP_r;

for i = 1:nM
  fprintf('[Evaluate] Model %d/%d: %s\n', i, nM, S.names{i});
  x = gdif_cache{i};
  Mv   = S.eval_mask;

  prd  = x(Mv) > cfg.decision_margin;
  yr   = S.y_raw_vec(Mv);
  ya   = S.y_adj_vec(Mv);

  m_raw = compute_metrics(prd, yr);
  m_adj = compute_metrics(prd, ya);

  % ---- MATH CHECKS: metrics identities ----
  TP_r(i)=m_raw.TP; FP_r(i)=m_raw.FP; TN_r(i)=m_raw.TN; FN_r(i)=m_raw.FN;
  TP_a(i)=m_adj.TP; FP_a(i)=m_adj.FP; TN_a(i)=m_adj.TN; FN_a(i)=m_adj.FN;
  ACC_r(i)=m_raw.ACC; ACC_a(i)=m_adj.ACC; PR_r(i)=m_raw.PREC; PR_a(i)=m_adj.PREC;
  RC_r(i)=m_raw.REC;  RC_a(i)=m_adj.REC;  F1_r(i)=m_raw.F1;   F1_a(i)=m_adj.F1;
  TPFP_r(i)=m_raw.TP_FP; TNFN_r(i)=m_raw.TN_FN; TPFP_a(i)=m_adj.TP_FP; TNFN_a(i)=m_adj.TN_FN;
end

S.results_table = table( ...
S.titles(:), nnz(S.eval_mask)*ones(nM,1), ...
TP_r, FP_r, TN_r, FN_r, TP_a, FP_a, TN_a, FN_a, ...
ACC_r, ACC_a, PR_r, PR_a, RC_r, RC_a, F1_r, F1_a, ...
TPFP_r, TNFN_r, TPFP_a, TNFN_a, ...
'VariableNames', {'Model','N_roi', ...
'TP_raw','FP_raw','TN_raw','FN_raw','TP_adj','FP_adj','TN_adj','FN_adj', ...
'ACC_raw','ACC_adj','PREC_raw','PREC_adj','REC_raw','REC_adj','F1_raw','F1_adj', ...
'TP_FP_raw','TN_FN_raw','TP_FP_adj','TN_FN_adj'});

%% ---- EXPECTED METRICS + (PIXEL/COMPONENT) BOOTSTRAP — REWRITTEN ----
fprintf('Analytic expected metrics with %s bootstrap...\n', cfg.uncertainty.bootstrap_mode);

nM = numel(S.names);
stat_list = {'SPEC','REC','PREC','ACC','F1'};

% ROI-vector setup
M = true(numel(idx),1);                % all ROI pixels participate equally
NM = nnz(M);
y_raw_M = double(S.y_raw_full(idx));   % already filled earlier
y_adj_M = double(S.y_adj_full(idx));

% sigma(Gmin) in ROI (or zeros if disabled)
if ~isempty(S.sigma_gmin), sigG_M = S.sigma_gmin(idx); else, sigG_M = zeros(NM,1,'like',y_raw_M); end

% Containers
EXP_raw = struct(); EXP_adj = struct();
CI_raw  = struct(); CI_adj  = struct();
CIi_raw = struct(); CIi_adj = struct();
for s = stat_list
    EXP_raw.(s{1}) = nan(nM,1);
    EXP_adj.(s{1}) = nan(nM,1);
    CI_raw.(s{1})  = nan(nM,2);
    CI_adj.(s{1})  = nan(nM,2);
    CIi_raw.(s{1}) = nan(nM,2);
    CIi_adj.(s{1}) = nan(nM,2);
end

alpha_main  = 1 - cfg.uncertainty.band_main;   % e.g., 0.10 for 90%
alpha_inner = 1 - cfg.uncertainty.band_inner;  % e.g., 0.50 for 50%

% RNG seed (if provided)
if isfield(cfg,'uncertainty') && isfield(cfg.uncertainty,'seed') && ~isempty(cfg.uncertainty.seed), rng(cfg.uncertainty.seed); end

% Use component bootstrap only if sinks are on and groups exist
use_comp_boot = strcmpi(cfg.uncertainty.bootstrap_mode,'component') && S.boot_groups.nC>0;

G_M = {};  % Pre-load groups in M-space (already built earlier as indices 1..NM)
if use_comp_boot, G_M = S.boot_groups.groups; end

for i = 1:nM   
    name_i = S.names{i};
    fprintf('Per-pixel σ\n');
    
    % --------- GΔ and σ(model) in M-space ----------
    gdif_M_all = gdif_cache{i}(M);

    % Per-pixel σ(model) for this model, aligned to M-space
    sigM_M = zeros(NM,1,'like',gdif_M_all);
    switch name_i
        case {'Hazzard','Martos','Shen','Stal'}
            if isfield(S.unc, name_i) && ~isempty(S.unc.(name_i))
                sigM_roi = S.unc.(name_i)(idx);    % ROI-vector
                sigM_M   = sigM_roi(M);            % M-space
            else
                sigM_M(:) = get_field_default(fallback_err, name_i, 0);
            end
        case 'Losing'
            if isfield(S,'bounds') && isfield(S.bounds,'Losing_min') && isfield(S.bounds,'Losing_max') ...
               && ~isempty(S.bounds.Losing_min) && ~isempty(S.bounds.Losing_max)
                half_roi = 0.5*abs(S.bounds.Losing_max(idx) - S.bounds.Losing_min(idx));
                sigM_M   = (half_roi(M)) / 1.96;   % σ from 95% range
            else
                sigM_M(:) = get_field_default(fallback_err, name_i, 0);
            end
        otherwise
            sigM_M(:) = get_field_default(fallback_err, name_i, 0);
    end

    % Pr[wet] in M-space, then finite mask V (per-model)
    pw_all = wet_prob_with_margin(gdif_M_all, sigM_M, sigG_M, cfg.decision_margin);
    V      = isfinite(pw_all) & isfinite(y_raw_M) & isfinite(y_adj_M);

    % Compress to finite subset for expectations & pixel bootstrap
    pw = pw_all(V);
    yr = y_raw_M(V);
    ya = y_adj_M(V);
    Nf = numel(pw);

    % --- Expected (non-bootstrap) confusion counts & metrics (RAW/ADJ) ---
    % RAW:
    ETP  = sum(pw .* yr);
    S1   = sum(pw);          % expected positives
    Ppos = sum(yr);          % actual positives
    EFP  = S1 - ETP;
    EFN  = Ppos - ETP;
    ETN  = (Nf - Ppos) - EFP;
    
    clamp = @(x) max(x, 0);
    [ETP,EFP,ETN,EFN] = deal(clamp(ETP), clamp(EFP), clamp(ETN), clamp(EFN));

    mR   = metrics_from_counts(ETP,EFP,ETN,EFN);
    EXP_raw.SPEC(i)=mR.SPEC; EXP_raw.REC(i)=mR.REC; EXP_raw.PREC(i)=mR.PREC;
    EXP_raw.ACC(i) =mR.ACC;  EXP_raw.F1(i) =mR.F1;

    % ADJ:
    ETPa = sum(pw .* ya);
    S1a  = S1;               % same S1 (depends only on pw)
    Pposa= sum(ya);
    EFPa = S1a - ETPa;
    EFNa = Pposa - ETPa;
    ETNa = (Nf - Pposa) - EFPa;
    mA   = metrics_from_counts(ETPa,EFPa,ETNa,EFNa);
    EXP_adj.SPEC(i)=mA.SPEC; EXP_adj.REC(i)=mA.REC; EXP_adj.PREC(i)=mA.PREC;
    EXP_adj.ACC(i) =mA.ACC;  EXP_adj.F1(i) =mA.F1;

    if startsWith(name_i,'Flat_'), continue; end

    % ---- Bootstrap (pixel or component) ----
    B  = cfg.uncertainty.n_boot;
    BR = struct(); BA = struct();
    for s = stat_list
        BR.(s{1}) = zeros(B,1);
        BA.(s{1}) = zeros(B,1);
    end

    if strcmpi(cfg.uncertainty.bootstrap_mode,'pixel')
        % ------- Pixel bootstrap over the finite subset (size = Nf) -------
        for b = 1:B
            rb = randi(Nf, Nf, 1);   % resample pixels with replacement
            % RAW
            ETPb = sum(pw(rb) .* yr(rb));
            S1b  = sum(pw(rb));
            Pposb= sum(yr(rb));
            EFPb = S1b - ETPb;
            EFNb = Pposb - ETPb;
            ETNb = (Nf - Pposb) - EFPb;
            mbR  = metrics_from_counts(ETPb,EFPb,ETNb,EFNb);
            % ADJ
            ETPab = sum(pw(rb) .* ya(rb));
            S1ab  = S1b;                 % same S1 under pixel sampling
            Pposab= sum(ya(rb));
            EFPab = S1ab - ETPab;
            EFNab = Pposab - ETPab;
            ETNab = (Nf - Pposab) - EFPab;
            mbA   = metrics_from_counts(ETPab,EFPab,ETNab,EFNab);
            % collect
            BR.SPEC(b)=mbR.SPEC; BR.REC(b)=mbR.REC; BR.PREC(b)=mbR.PREC; BR.ACC(b)=mbR.ACC; BR.F1(b)=mbR.F1;
            BA.SPEC(b)=mbA.SPEC; BA.REC(b)=mbA.REC; BA.PREC(b)=mbA.PREC; BA.ACC(b)=mbA.ACC; BA.F1(b)=mbA.F1;
        end

        fprintf('[boot pixel] %-10s  N=%d  B=%d  mean(pw)=%.3f  sd(pw)=%.3f  frac~{0,1}=%.2f\n', ...
                name_i, Nf, B, mean(pw,'omitnan'), std(pw,'omitnan'), ...
                mean((pw<1e-3)|(pw>1-1e-3),'omitnan'));
    else
        % ------- Component bootstrap (resample connected sinks) -------
        if isempty(G_M), continue; end

        % Filter each component's member indices by the per-model finite mask V
        GV = cellfun(@(g) g(V(g)), G_M, 'UniformOutput', false);
        GV = GV(~cellfun('isempty', GV));
        nC = numel(GV);
        if nC == 0, continue; end

        for b = 1:B
            mC = max(1, round(cfg.uncertainty.mfrac * nC));
            rbC = randi(nC, mC, 1);          % sample mC components with replacement            
            pick = vertcat(GV{rbC});         % concatenated member pixels (finite by construction)

            % RAW
            ETPb = sum(pw_all(pick) .* y_raw_M(pick));
            S1b  = sum(pw_all(pick));
            Pposb= sum(y_raw_M(pick));
            EFPb = S1b - ETPb;
            EFNb = Pposb - ETPb;
            ETNb = (numel(pick) - Pposb) - EFPb;
            mbR  = metrics_from_counts(ETPb,EFPb,ETNb,EFNb);

            % ADJ
            ETPab = sum(pw_all(pick) .* y_adj_M(pick));
            S1ab  = S1b;
            Pposab= sum(y_adj_M(pick));
            EFPab = S1ab - ETPab;
            EFNab = Pposab - ETPab;
            ETNab = (numel(pick) - Pposab) - EFPab;
            mbA   = metrics_from_counts(ETPab,EFPab,ETNab,EFNab);

            BR.SPEC(b)=mbR.SPEC; BR.REC(b)=mbR.REC; BR.PREC(b)=mbR.PREC; BR.ACC(b)=mbR.ACC; BR.F1(b)=mbR.F1;
            BA.SPEC(b)=mbA.SPEC; BA.REC(b)=mbA.REC; BA.PREC(b)=mbA.PREC; BA.ACC(b)=mbA.ACC; BA.F1(b)=mbA.F1;
        end

        fprintf('[boot comp ] %-10s  nComp=%d  B=%d  mean(pw)=%.3f  sd(pw)=%.3f  frac~{0,1}=%.2f\n', ...
            name_i, nC, B, mean(pw_all(V),'omitnan'), std(pw_all(V),'omitnan'), ...
            mean((pw_all(V)<1e-3)|(pw_all(V)>1-1e-3),'omitnan'));
    end

    % ---- Quantiles → outer (band_main) and inner (band_inner) CIs ----
    for s = stat_list
        vR = BR.(s{1}); vA = BA.(s{1});
        if ~isempty(vR)
            CI_raw.(s{1})(i,:)  = prctile(vR, 100*[alpha_main/2, 1 - alpha_main/2]);
            CIi_raw.(s{1})(i,:) = prctile(vR, 100*[alpha_inner/2, 1 - alpha_inner/2]);
        end
        if ~isempty(vA)
            CI_adj.(s{1})(i,:)  = prctile(vA, 100*[alpha_main/2, 1 - alpha_main/2]);
            CIi_adj.(s{1})(i,:) = prctile(vA, 100*[alpha_inner/2, 1 - alpha_inner/2]);
        end
    end
end

% ---- Pack results ----
S.Gstats = struct();
for s = stat_list
    S.Gstats.([s{1} '_raw']) = struct('EXP', EXP_raw.(s{1}), 'CI', CI_raw.(s{1}), 'CIi', CIi_raw.(s{1}));
    S.Gstats.([s{1} '_adj']) = struct('EXP', EXP_adj.(s{1}), 'CI', CI_adj.(s{1}), 'CIi', CIi_adj.(s{1}));
end

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

%% -------------------- COLORED-CI PLOTS (single & generalized) ------------
assert(isfield(S,'Gstats'), 'Gstats missing — run uncertainty/bootstraps first');
try
    px = upper('SPEC'); py = upper('REC');      % Basic field presence check for this pair
    needX = [px '_adj']; needY = [py '_adj'];
   
    GX = S.Gstats.(needX); GY = S.Gstats.(needY);  % Extra sanity: CI order & finite EXP (vectorized quick check)
    
    assert(all(isfinite(GX.EXP) | isnan(GX.EXP)), '%s ADJ EXP has non-finite non-NaN entries', px);
    assert(all(isfinite(GY.EXP) | isnan(GY.EXP)), '%s ADJ EXP has non-finite non-NaN entries', py);

    % ---- Plot panel (ADJ) ----
    hf = plot_colored_ci_pair(S, cfg, 'SPEC', 'REC', 'adj');
    set(findall(hf, '-property', 'FontSize'), 'FontSize', cfg.font_size);
    clim([35 60]);
    out = fullfile(cfg.outdir, sprintf('ci2_%s_%s_ADJ%s.png', ...
        upper('SPEC'), upper('REC'), region_suffix(cfg)));
    exportgraphics(hf, out, 'Resolution', 300);

catch ME
    warning('2D plot %s-%s ADJ failed: %s', 'SPEC', 'REC', ME.message);
end

%% ===== READOUTS (CI sanity + quick table) ===========================
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

%% -------------------- MAPS -------------------------------------
fprintf('Rendering maps (GΔ masked by slow flow; sinks overlaid as points)...\n');

% ---- Masks & styles ----
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
%L_ref = L;  % ensure defined even if try/catch path varies
fprintf('Diverging limits set to [-%.3g, %.3g]\n', L, L);

% ---- GΔ maps for each non-flat model showing predictions ----
nM = numel(S.names);
for ii = 1:nM
  name = S.names{ii};
  title_i = S.titles{ii};
  if startsWith(name,'Flat_'), continue; end

  zd_path = fullfile(cfg.cache_dir, sprintf('Gdif_full_%s.mat', name));
  if isfile(zd_path)
      Ld = load(zd_path,'Gdif_full');  Zd = Ld.Gdif_full;
  else
      [~, ~, Zd] = processGHF(S.models.(name), S.Ts, S.Mb, S.H);
      try Gdif_full = Zd; save(zd_path,'Gdif_full','-v7.3'); catch, end
  end

  % --- after you have Zd and title_i ---
  Zroi = Zd(S.roi_mask & isfinite(Zd));
  if isempty(Zroi), Zroi = Zd(isfinite(Zd)); end
  L = max(abs(prctile(double(Zroi),[2 98]))); if ~isfinite(L) || L<=0, L = 50; end
    
  fig = figure('Visible','off','Color','w','Position',[100 100 900 800]);
    
  % Top/only axis: GΔ overlay
  axZ = axes('Parent',fig,'Color','none'); hold(axZ,'on'); axis(axZ,'image'); box(axZ,'on');
  set(axZ,'YDir','normal','Layer','top');

  him = imagesc(axZ, S.X_km(1,:), S.Y_km(:,1), Zd);
  set(him, 'AlphaData', cfg.overlay_alpha * (S.viz_mask & isfinite(Zd)), ...
             'AlphaDataMapping', 'none');
  colormap(axZ, diverging_cmap(256));
  clim(axZ, [-L L]);
  cb = colorbar(axZ); ylabel(cb, 'G_\Delta (mW m^{-2})');
    
  xlim([0 700]);
  ylim([-200 350]);
    
  try
      step = 500;
      Hmin = floor(min(S.H(:),[],'omitnan')/step)*step;
      Hmax = ceil( max(S.H(:),[],'omitnan')/step)*step;
      if isfinite(Hmin) && isfinite(Hmax) && Hmax>Hmin
          levels = Hmin:step:Hmax;
          contour(axZ, S.X_km, S.Y_km, S.H, levels, 'LineColor',[0.3 0.3 0.3], 'LineWidth',0.5);
      end
  catch
  end
    
    mask_sinks = S.roi_mask & S.sink_mask & isfinite(Zd);
    
    Mw = mask_sinks & (Zd > 0);
    Md = mask_sinks & (Zd <= 0);
    [wy, wx] = find(Mw);
    [dy, dx] = find(Md);
    scatter(axZ, S.X_km(1,wx), S.Y_km(wy,1), 6, '.', 'MarkerEdgeColor',[0.83 0.23 0.23], 'MarkerEdgeAlpha',0.95);
    scatter(axZ, S.X_km(1,dx), S.Y_km(dy,1), 6, '.', 'MarkerEdgeColor',[0.20 0.40 0.90], 'MarkerEdgeAlpha',0.95);
    
    title(axZ, sprintf('G_\\Delta (diverging, ROI-masked) + sinks — %s', title_i));
    xlabel(axZ,'X (km)'); ylabel(axZ,'Y (km)');
    set(findall(fig, '-property', 'FontSize'), 'FontSize', cfg.font_size);
    
    outbase = sprintf('map_Gdif_%s%s.png', safe(name), region_suffix(cfg));
    outfile = fullfile(cfg.outdir, outbase);
    exportgraphics(fig, outfile, 'Resolution', 300); 
    clear Zd Ld
end

% ---- GΔ maps for each non-flat model showing where predictions are wrong ----
nM = numel(S.names);
for ii = 1:nM
  name = S.names{ii};
  title_i = S.titles{ii};
  if startsWith(name,'Flat_'), continue; end

  zd_path = fullfile(cfg.cache_dir, sprintf('Gdif_full_%s.mat', name));
  if isfile(zd_path)
      Ld = load(zd_path,'Gdif_full');  Zd = Ld.Gdif_full;
  else
      [~, ~, Zd] = processGHF(S.models.(name), S.Ts, S.Mb, S.H);
      try Gdif_full = Zd; save(zd_path,'Gdif_full','-v7.3'); catch, end
  end

  % --- after you have Zd and title_i ---
  Zroi = Zd(S.roi_mask & isfinite(Zd));
  if isempty(Zroi), Zroi = Zd(isfinite(Zd)); end
  L = max(abs(prctile(double(Zroi),[2 98]))); if ~isfinite(L) || L<=0, L = 50; end
    
  fig = figure('Visible','off','Color','w','Position',[100 100 900 800]);
    
  % Top/only axis: GΔ overlay
  axZ = axes('Parent',fig,'Color','none'); hold(axZ,'on'); axis(axZ,'image'); box(axZ,'on');
  set(axZ,'YDir','normal','Layer','top');

  him = imagesc(axZ, S.X_km(1,:), S.Y_km(:,1), Zd);
  set(him, 'AlphaData', cfg.overlay_alpha * (S.viz_mask & isfinite(Zd)), ...
             'AlphaDataMapping', 'none');
  colormap(axZ, diverging_cmap(256));
  clim(axZ, [-L L]);
  cb = colorbar(axZ); ylabel(cb, 'G_\Delta (mW m^{-2})');
    
  xlim([0 700]);
  ylim([-200 350]);
    
  try
      step = 500;
      Hmin = floor(min(S.H(:),[],'omitnan')/step)*step;
      Hmax = ceil( max(S.H(:),[],'omitnan')/step)*step;
      if isfinite(Hmin) && isfinite(Hmax) && Hmax>Hmin
          levels = Hmin:step:Hmax;
          contour(axZ, S.X_km, S.Y_km, S.H, levels, 'LineColor',[0.3 0.3 0.3], 'LineWidth',0.5);
      end
  catch
  end
    
  % ---- Sinks colored by prediction correctness (ADJ labels) ----
  % ---- Sinks colored by correctness (TP/TN/FP/FN) ----
  pred_full  = (Zd > cfg.decision_margin);   % model prediction at grid cells
  truth_full = S.y_adj_full;                 % or S.y_raw_full for raw labels

  mask_sinks = S.roi_mask & S.sink_mask & isfinite(Zd);

  TP = mask_sinks &  pred_full &  truth_full;   TN = mask_sinks & ~pred_full & ~truth_full;
  FP = mask_sinks &  pred_full & ~truth_full;   FN = mask_sinks & ~pred_full &  truth_full;

  [iy,ix] = find(TP);
  hTP = scatter(axZ, S.X_km(1,ix), S.Y_km(iy,1), 8, '.', ...
      'MarkerEdgeColor',[0.10 0.70 0.10], 'MarkerEdgeAlpha',0.95);

  [iy,ix] = find(TN);
  hTN = scatter(axZ, S.X_km(1,ix), S.Y_km(iy,1), 8, '.', ...
      'MarkerEdgeColor',[0.20 0.40 0.90], 'MarkerEdgeAlpha',0.95);

  [iy,ix] = find(FP);
  hFP = scatter(axZ, S.X_km(1,ix), S.Y_km(iy,1), 8, '.', ...
      'MarkerEdgeColor',[0.90 0.55 0.10], 'MarkerEdgeAlpha',0.95);

  [iy,ix] = find(FN);
  hFN = scatter(axZ, S.X_km(1,ix), S.Y_km(iy,1), 8, '.', ...
      'MarkerEdgeColor',[0.85 0.20 0.60], 'MarkerEdgeAlpha',0.95);

  legend(axZ, [hTP hTN hFP hFN], {'TP','TN','FP','FN'}, 'Location','southoutside');

  xlabel(axZ,'X (km)'); ylabel(axZ,'Y (km)');
  set(findall(fig, '-property', 'FontSize'), 'FontSize', cfg.font_size);
    
    % build a name and guarantee an extension
  outbase = sprintf('map_Gdif_correct_%s%s.png', safe(name), region_suffix(cfg));
  outfile = fullfile(cfg.outdir, outbase);
  try
      exportgraphics(fig, outfile, 'Resolution', 300);
  catch
      print(fig, outfile, '-dpng', '-r300');
  end

    clear Zd Ld
end

% ---- Gmin map with sinks colored by specularity (red=high, blue=low) ----
[loGmin, hiGmin] = robust_clim_fast(Gmin_full, S.viz_mask);   % 2–98% over viz
fig = figure('Visible','off','Color','w','Position',[100 100 900 800]);

% Background: H grayscale, overlay: Gmin sequential
axH = axes('Parent',fig); hold(axH,'on'); axis(axH,'image'); box(axH,'on');
set(axH,'YDir','normal','Layer','bottom');
Hfinite = S.H(isfinite(S.H)); Hlo = prctile(double(Hfinite),2); Hhi = prctile(double(Hfinite),98);
imagesc(axH, S.X_km(1,:), S.Y_km(:,1), mat2gray(double(S.H),[Hlo Hhi])); colormap(axH, gray);
set(axH,'XTick',[],'YTick',[],'XColor','none','YColor','none','Box','off');

axZ = axes('Parent',fig,'Color','none'); hold(axZ,'on'); axis(axZ,'image'); box(axZ,'on');
set(axZ,'YDir','normal','Layer','top');
hG = imagesc(axZ, S.X_km(1,:), S.Y_km(:,1), Gmin_full);
set(hG,'AlphaData', double(S.viz_mask) * cfg.overlay_alpha, 'AlphaDataMapping','none');
colormap(axZ, parula(256)); clim(axZ, [loGmin hiGmin]);
cb = colorbar(axZ); ylabel(cb, 'G_{min} (mW m^{-2})');

% Axis frame and extents
xlim(axZ,[0 700]); ylim(axZ,[-200 350]);
xlabel(axZ,'X (km)'); ylabel(axZ,'Y (km)');
setappdata(fig,'linkXY',linkprop([axH,axZ],{'XLim','YLim'})); axis(axZ,'image');

%title(axZ, sprintf('G_{min} with sinks (Q \\ge %.2f = red, < %.2f = blue)', Qthr, Qthr));
set(findall(fig,'-property','FontSize'),'FontSize',cfg.font_size);

outfile = fullfile(cfg.outdir, sprintf('map_Gmin%s.png', region_suffix(cfg)));
   exportgraphics(fig, outfile, 'Resolution', 300);  

end

%% ===================== Helpers =====================
function s = region_suffix(cfg)
    s = '';
    try
        if isfield(cfg,'region') && isfield(cfg.region,'mode') && ~isempty(cfg.region.mode)
            s = ['_REG-' upper(char(cfg.region.mode))];
        end
    catch
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
  if nargin < 2 || isempty(M), Zs = Z(isfinite(Z));
      else, Zs = Z(isfinite(Z) & M);
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

% MD5 of a numeric array (size + data), robust cache key
function h = md5_of_array(A)
    if isempty(A), h = repmat('0',1,32); return; end
    bs   = getByteStreamFromArray(A);
    md   = java.security.MessageDigest.getInstance('MD5');
    md.update(uint8(bs));
    raw  = typecast(md.digest(),'uint8');
    h    = lower(reshape(dec2hex(raw,2).',1,[]));
end

% Finite-difference uncertainty propagation for Gmin
% sigma(Gmin)^2 ≈ (∂G/∂Ts)^2 σ_Ts^2 + (∂G/∂Mb)^2 σ_Mb^2 + (∂G/∂H)^2 σ_H^2
function [sigma_map, sigma_scalar] = estimate_sigma_gmin_fd(S, cfg)
  
    key = sprintf('%s_%s_%s', md5_of_array(S.H), md5_of_array(S.Ts), md5_of_array(S.Mb));
    cache_path = fullfile(cfg.cache_dir, ['sigma_gmin_fd_' key '.mat']);
    if isfile(cache_path)
        L = load(cache_path, 'sigma_map', 'sigma_scalar');
        sigma_map = L.sigma_map; sigma_scalar = L.sigma_scalar; return;
    end

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

  if isfield(S,'roi_mask') && ~isempty(S.roi_mask), vv = sigma_map(S.roi_mask & isfinite(sigma_map));
      else, vv = sigma_map(isfinite(sigma_map));
  end
  sigma_scalar = sqrt(mean(double(vv).^2, 'omitnan'));

  try save(cache_path, 'sigma_map','sigma_scalar','-v7.3'); catch; end
end

function val = get_field_default(s, fld, defaultVal)
    if isstruct(s) && isfield(s, fld) && ~isempty(s.(fld)), val = s.(fld);
        else, val = defaultVal;
    end
end

function m = metrics_from_counts(TP,FP,TN,FN)
    N = TP + TN + FP + FN + eps;
    PREC = TP / max(TP + FP, eps);
    SPEC = TN / max(TN + FP, eps);
    REC  = TP / max(TP + FN, eps);
    ACC  = (TP + TN) / N;
    F1   = 2 * PREC * REC / max(PREC + REC, eps);
    m = struct('PREC',PREC,'SPEC',SPEC,'REC',REC,'ACC',ACC,'F1',F1);
end

function lim = stat_axis_limits(stat_name)
    s = upper(stat_name);
    switch s
        case 'MCC'
            lim = [-1 1];
        otherwise
            lim = [0 1];
    end
end

function f = plot_colored_ci_pair(S, cfg, statX, statY, mode)
% plot colored CI rectangles for (statX, statY)

    statX = upper(statX); statY = upper(statY); mode = lower(mode);
    GX = S.Gstats.([statX '_' mode]);
    GY = S.Gstats.([statY '_' mode]);

    X  = GX.EXP;   Xlohi  = GX.CI;   Xlohi_i = GX.CIi;
    Y  = GY.EXP;   Ylohi  = GY.CI;   Ylohi_i = GY.CIi;

    nM = numel(S.names);
    is_flat = startsWith(S.names,'Flat_');
    colors = lines(max(1,nM));
    epsw = 1e-6;

    f = figure('Color','w','Position',[120 120 900 820]);
    ax = axes(f); hold(ax,'on'); grid(ax,'on'); box(ax,'on'); ax.Layer='top';
    xlabel(ax, sprintf('Fraction of "dry" sinks correctly predicted'));
    ylabel(ax, sprintf('Fraction of "wet" sinks correctly predicted'));
    xlim(ax, stat_axis_limits(statX));
    ylim(ax, stat_axis_limits(statY));

    If = find(is_flat);        %Flat GHF sweep (gradient by flat value in mW m^-2)
    xs = X(If); ys = Y(If); vals = S.model_cvals(If); 
    good = isfinite(xs) & isfinite(ys) & isfinite(vals);
    xs = xs(good); ys = ys(good); vals = vals(good);

    if numel(vals) >= 2
        [vals_s, ord] = sort(vals);
        xs_s = xs(ord); ys_s = ys(ord);

        if isfield(cfg.flat,'interp_enable') && cfg.flat.interp_enable && numel(vals_s) >= 2
            npts = max(20, cfg.flat.interp_points);
            [u, iu] = unique(vals_s);                  % ensure strictly increasing param
            vi = linspace(min(u), max(u), npts);       % parameter = flat value (mW m^-2)
            mth = 'spline';
            if isfield(cfg.flat,'interp_method') && ~isempty(cfg.flat.interp_method), mth = cfg.flat.interp_method; end
            xi = interp1(u, xs_s(iu), vi, mth, 'extrap');
            yi = interp1(u, ys_s(iu), vi, mth, 'extrap');

            plot_gradient_curve(ax, xi, yi, vi, 20);
        else
            plot_gradient_curve(ax, xs_s, ys_s, vals_s, 20);
        end

        colormap(ax, parula(256));
        clim(ax, [30 70]);
        grid off;
        cb = colorbar(ax);
        ylabel(cb, 'Flat GHF (mW m^{-2})');
    end

    for i = 1:nM    % Non-flat models: colored CI envelopes + markers + labels
        if is_flat(i), continue; end
        if ~all(isfinite([X(i) Y(i)])), continue; end

        % outer band
        if all(isfinite([Xlohi(i,:) Ylohi(i,:)]))
            x1 = Xlohi(i,1); x2 = Xlohi(i,2); if x2<=x1, x2 = x1+epsw; end
            y1 = Ylohi(i,1); y2 = Ylohi(i,2); if y2<=y1, y2 = y1+epsw; end
            patch(ax, [x1 x2 x2 x1], [y1 y1 y2 y2], colors(i,:), ...
                  'FaceAlpha',0.3, 'EdgeColor','none', 'HitTest','off');
        end

        % inner band
        if all(isfinite([Xlohi_i(i,:) Ylohi_i(i,:)]))
            x1 = Xlohi_i(i,1); x2 = Xlohi_i(i,2); if x2<=x1, x2 = x1+epsw; end
            y1 = Ylohi_i(i,1); y2 = Ylohi_i(i,2); if y2<=y1, y2 = y1+epsw; end
            patch(ax, [x1 x2 x2 x1], [y1 y1 y2 y2], colors(i,:), ...
                  'FaceAlpha',0.6, 'EdgeColor','none', 'HitTest','off');
        end

        plot(ax, X(i), Y(i), 'o', 'MarkerFaceColor',colors(i,:), ...
             'MarkerEdgeColor','k', 'MarkerSize',6, 'HandleVisibility','off');
        if i <= numel(S.titles)
            dx = 0.01 * range(xlim(ax)); dy = 0.01 * range(ylim(ax));
            text(ax, X(i)+dx, Y(i)+dy, S.titles{i}, ...
                 'FontSize',20, 'Color',colors(i,:), 'Interpreter','none');
        end
    end
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

% Pr[(M - Gmin) > m] when (M - Gmin) ~ N(mu, sig_m^2 + sig_g^2)
function pwet = wet_prob_with_margin(gdif_vec, sig_m, sig_g, margin)
  if isempty(sig_g), sig_g = 0; end
  if nargin < 4 || isempty(margin), margin = 0; end
  s2 = sig_m.^2 + sig_g.^2;
  pwet = NaN(size(gdif_vec));
  Z = (s2==0);
  % exact zero-σ: hard step with tie=0.5 at gdif==margin
  pwet(Z) = double(gdif_vec(Z) > margin) + 0.5*double(gdif_vec(Z) == margin);
  NZ = ~Z;
  if any(NZ(:))
      s = sqrt(max(s2(NZ), eps));
      pwet(NZ) = 0.5*(1 + erf(((gdif_vec(NZ) - margin)./s)/sqrt(2)));
  end
end

% Split grid into ISPB/OSPB using a vertical scanline against the boundary.
function [mask_ISPB, mask_OSPB, xb_row] = build_spb_split_masks(X_km, Y_km, Bxy_km)
    yrows = Y_km(:,1);                % row y (km)
    xb_row = nan(numel(yrows),1);     % boundary x at each row y

    bx = Bxy_km(:,1); by = Bxy_km(:,2);
    good = isfinite(bx) & isfinite(by);
    bx = bx(good); by = by(good);

    for j = 1:numel(by)-1      % Accumulate per-row intersections; if multiple, use the max-x (outer edge)
        x1 = bx(j);  y1 = by(j);
        x2 = bx(j+1);y2 = by(j+1);
        if ~isfinite(x1) || ~isfinite(y1) || ~isfinite(x2) || ~isfinite(y2), continue; end

        lo = min(y1,y2); hi = max(y1,y2);        % rows that intersect this segment's y-span
        I  = find(yrows >= lo & yrows <= hi);
        if isempty(I), continue; end

        if y2 ~= y1
            t = (yrows(I) - y1) ./ (y2 - y1);
            xint = x1 + t .* (x2 - x1);
            if any(isfinite(xint))            % keep the rightmost (max x) intersection per row
                if any(isfinite(xb_row(I))), xb_row(I) = max(xb_row(I), xint);
                    else, xb_row(I) = xint;
                end
            end
        else
            % Horizontal segment; snap the exact row(s) if any
            I0 = find(abs(yrows - y1) <= eps(max(abs(yrows))));
            if ~isempty(I0)
                xr = max(x1,x2);
                if any(isfinite(xb_row(I0))), xb_row(I0) = max(xb_row(I0), xr);
                    else, xb_row(I0) = xr;
                end
            end
        end
    end

    xb_row = fillmissing(xb_row,'nearest');      % Fill gaps by nearest row; clamp to grid x-range
    x_min = min(X_km(1,:)); x_max = max(X_km(1,:));
    xb_row = min(max(xb_row, x_min), x_max);
 
    mask_ISPB = bsxfun(@le, X_km, xb_row);  % Broadcast rows to full grid and build masks
    mask_OSPB = ~mask_ISPB;
end

function M = build_rect_mask(X_km, Y_km, rectCfg)
    M = true(size(X_km));                 % default = all pass
    if nargin<3 || isempty(rectCfg) || ~isstruct(rectCfg), return; end
    if ~isfield(rectCfg,'enable') || ~rectCfg.enable, return; end

    rects = [];
    if isfield(rectCfg,'rects') && ~isempty(rectCfg.rects)
        rects = rectCfg.rects;
    elseif isfield(rectCfg,'xlim') && ~isempty(rectCfg.xlim) && ...
           isfield(rectCfg,'ylim') && ~isempty(rectCfg.ylim)
        rects = [rectCfg.xlim(:).' rectCfg.ylim(:).'];
    end
   
    if isempty(rects), return; end
    if size(rects,2) ~= 4, error('cfg.rect.rects must be N×4 rows: [xmin xmax ymin ymax]'); end

    keep = false(size(X_km));
    for r = 1:size(rects,1)
        xmin = min(rects(r,1), rects(r,2));
        xmax = max(rects(r,1), rects(r,2));
        ymin = min(rects(r,3), rects(r,4));
        ymax = max(rects(r,3), rects(r,4));
        keep = keep | (X_km >= xmin & X_km <= xmax & Y_km >= ymin & Y_km <= ymax);
    end

    mode = 'include';
    if isfield(rectCfg,'mode') && ~isempty(rectCfg.mode), mode = lower(string(rectCfg.mode)); end
    if strcmp(mode,'exclude'), M = ~keep;
        else, M = keep;
    end
end