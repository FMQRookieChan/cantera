function results = ndodecane_combustor_2d(userConfig)
    %% Two-dimensional n-dodecane combustor model on a fixed 180 x 70 mesh
    %
    % Compare a single fuel inlet and a dual fuel inlet configuration on the
    % same structured mesh. Arrays are indexed as [radial, axial] = [70, 180].
    %
    % The V2500 / PW1100G inlet values are representative engineering
    % assumptions and are not certified manufacturer data.
    %
    % Chemistry architecture:
    %   - Preferred: Cantera MATLAB callback / lookup using
    %     data/nDodecane_Reitz.yaml (phase nDodecane_IG by default)
    %   - Fallback: bounded Arrhenius surrogate when Cantera toolbox is absent
    %
    % This script is a compact advection-diffusion-reaction demonstrator and
    % does not claim full CFD fidelity.
    %
    % Example:
    %   results = ndodecane_combustor_2d(struct('engine', 'PW1100G'));
    %
    % Requires: cantera >= 3.2.0 for Cantera chemistry mode
    %
    % .. tags:: Matlab, combustion, finite volume, plotting

    if nargin < 1
        userConfig = struct();
    end

    cfg = applyOverrides(defaultConfig(), userConfig);
    mesh = buildMesh(cfg);
    engine = getEngineCase(cfg.engine);

    chemistry = initializeChemistry(cfg);

    results = struct();
    results.metadata = struct(
        'meshShape', [mesh.nr, mesh.nx], ...
        'indexing', '[radial, axial]', ...
        'engine', engine.name, ...
        'engineAssumptions', engine.assumptionNote, ...
        'mechanismPath', chemistry.mechanismPath, ...
        'chemistryMode', chemistry.mode, ...
        'phase', cfg.phase);

    results.single = runCase(mesh, cfg, engine, chemistry, 'single');
    results.dual = runCase(mesh, cfg, engine, chemistry, 'dual');

    if cfg.makePlots
        plotCaseFields(results, mesh, engine, cfg);
    end
end

function cfg = defaultConfig()
    cfg = struct();

    % Fixed mesh requirement
    cfg.nx = 180;
    cfg.nr = 70;
    cfg.length = 0.45;   % m
    cfg.radius = 0.06;   % m

    cfg.engine = 'V2500';

    % Representative overall operating settings
    cfg.phi = 0.38;
    cfg.YO2In = 0.233;
    cfg.cp = 1150.0;     % J/(kg K)
    cfg.Rgas = 287.0;    % J/(kg K)
    cfg.LHV = 43.0e6;    % J/kg
    cfg.gamma = 1.33;

    % Transport / numerics
    cfg.DFuel = 2.0e-5;  % m^2/s
    cfg.DO2 = 2.5e-5;    % m^2/s
    cfg.DNO = 2.5e-5;    % m^2/s
    cfg.alpha = 2.0e-5;  % m^2/s
    cfg.cfl = 0.45;
    cfg.dtMax = 5.0e-5;
    cfg.tEnd = 0.012;

    % Numerical safeguards
    cfg.minTemperature = 250.0;
    cfg.maxTemperature = 3200.0;
    cfg.pressureFloor = 2.0e5;
    cfg.maxMassFraction = 0.995;
    cfg.maxFuelRate = 1.0e3;

    % Pressure model coefficients
    cfg.pressureLossK = 0.08;
    cfg.pressureHeatK = 0.015;

    % Mechanism / chemistry
    cfg.mechanism = 'nDodecane_Reitz.yaml';
    cfg.phase = 'nDodecane_IG';
    cfg.forceFallbackChemistry = false;

    % Single-inlet fuel boundary region (radial fraction)
    cfg.singleRegion = struct('rMinFrac', 0.30, 'rMaxFrac', 0.55);

    % Dual internal injection regions
    cfg.dualSplit = [0.55, 0.45];
    cfg.dualRegions = {
        struct('xMinFrac', 0.15, 'xMaxFrac', 0.28, 'rMinFrac', 0.22, 'rMaxFrac', 0.42), ...
        struct('xMinFrac', 0.56, 'xMaxFrac', 0.70, 'rMinFrac', 0.50, 'rMaxFrac', 0.72) ...
    };

    cfg.makePlots = true;
end

function cfg = applyOverrides(cfg, userConfig)
    fn = fieldnames(userConfig);
    for i = 1:numel(fn)
        cfg.(fn{i}) = userConfig.(fn{i});
    end

    if cfg.nx ~= 180 || cfg.nr ~= 70
        error('Mesh is fixed to 180 axial x 70 radial cells.');
    end

    if abs(sum(cfg.dualSplit) - 1.0) > 1e-12
        error('dualSplit must sum to one.');
    end

    if numel(cfg.dualRegions) ~= 2
        error('dualRegions must contain exactly two injection regions.');
    end
end

function mesh = buildMesh(cfg)
    mesh = struct();
    mesh.nx = cfg.nx;
    mesh.nr = cfg.nr;
    mesh.dx = cfg.length / cfg.nx;
    mesh.dr = cfg.radius / cfg.nr;

    xCenters = ((1:cfg.nx) - 0.5) * mesh.dx;
    rCenters = ((1:cfg.nr) - 0.5) * mesh.dr;
    [mesh.X, mesh.R] = meshgrid(xCenters, rCenters);

    % Axisymmetric ring areas at x-inlet face and cell volumes
    rInner = ((1:cfg.nr) - 1) * mesh.dr;
    rOuter = (1:cfg.nr) * mesh.dr;
    mesh.inletRingArea = pi * (rOuter.^2 - rInner.^2)';
    mesh.cellVolume = repmat(mesh.inletRingArea * mesh.dx, 1, cfg.nx);

    if ~isequal(size(mesh.X), [70, 180])
        error('Mesh array size must be [70, 180] using [radial, axial] indexing.');
    end
end

function engine = getEngineCase(name)
    switch upper(string(name))
        case "V2500"
            engine = struct(...
                'key', 'V2500', ...
                'name', 'IAE V2500 representative combustor inlet', ...
                'Tin', 700.0, ...
                'Pin', 14.5e5, ...
                'mdotAir', 45.0, ...
                'assumptionNote', 'Representative engineering assumptions (not OEM-certified).');

        case "PW1100G"
            engine = struct(...
                'key', 'PW1100G', ...
                'name', 'PW1100G representative combustor inlet', ...
                'Tin', 730.0, ...
                'Pin', 16.0e5, ...
                'mdotAir', 50.0, ...
                'assumptionNote', 'Representative engineering assumptions (not OEM-certified).');

        otherwise
            error('Unknown engine case. Use V2500 or PW1100G.');
    end
end

function chemistry = initializeChemistry(cfg)
    chemistry = struct();
    chemistry.mode = 'fallback-arrhenius';
    chemistry.useCantera = false;
    chemistry.hasMechanismNO = false;
    chemistry.mechanismPath = resolveMechanismPath(cfg.mechanism);
    chemistry.phase = cfg.phase;

    % Fallback model constants
    chemistry.fallbackA = 3.5e3;
    chemistry.fallbackE = 8.2e4;
    chemistry.fallbackStoichO2 = 3.48;

    if cfg.forceFallbackChemistry
        return;
    end

    if exist('ct.Solution', 'class') ~= 8
        return;
    end

    try
        gas = ct.Solution(chemistry.mechanismPath, cfg.phase, 'none');
        fuelName = findSpeciesCaseInsensitive(gas, 'c12h26');
        o2Name = findSpeciesCaseInsensitive(gas, 'o2');
        n2Name = findSpeciesCaseInsensitive(gas, 'n2');

        if isempty(fuelName) || isempty(o2Name) || isempty(n2Name)
            return;
        end

        chemistry.useCantera = true;
        chemistry.mode = 'cantera-lookup-callback';
        chemistry.gas = gas;
        chemistry.fuelName = fuelName;
        chemistry.o2Name = o2Name;
        chemistry.n2Name = n2Name;
        chemistry.noName = findSpeciesCaseInsensitive(gas, 'no');
        chemistry.hasMechanismNO = ~isempty(chemistry.noName);

        chemistry.fuelIndex = gas.speciesIndex(chemistry.fuelName);
        chemistry.o2Index = gas.speciesIndex(chemistry.o2Name);
        chemistry.n2Index = gas.speciesIndex(chemistry.n2Name);
        if chemistry.hasMechanismNO
            chemistry.noIndex = gas.speciesIndex(chemistry.noName);
        else
            chemistry.noIndex = -1;
        end

        mw = gas.molecularWeights;
        chemistry.MWFuel = mw(chemistry.fuelIndex);
        chemistry.MWO2 = mw(chemistry.o2Index);
        chemistry.MWN2 = mw(chemistry.n2Index);
        if chemistry.hasMechanismNO
            chemistry.MWNO = mw(chemistry.noIndex);
        else
            chemistry.MWNO = 30.006;
        end

        chemistry.fallbackStoichO2 = (17.75 * chemistry.MWO2) / chemistry.MWFuel;

    catch
        chemistry.useCantera = false;
        chemistry.mode = 'fallback-arrhenius';
    end
end

function pathOut = resolveMechanismPath(pathIn)
    if isstring(pathIn)
        pathIn = char(pathIn);
    end

    if isfile(pathIn)
        pathOut = pathIn;
        return;
    end

    thisFile = mfilename('fullpath');
    repoRoot = fileparts(fileparts(fileparts(thisFile)));
    candidates = {
        fullfile(pwd, pathIn), ...
        fullfile(fileparts(thisFile), pathIn), ...
        fullfile(repoRoot, pathIn), ...
        fullfile(repoRoot, 'data', pathIn) ...
    };

    for i = 1:numel(candidates)
        if isfile(candidates{i})
            pathOut = candidates{i};
            return;
        end
    end

    error('Unable to resolve mechanism path: %s', pathIn);
end

function name = findSpeciesCaseInsensitive(gas, species)
    names = gas.speciesNames;
    idx = find(strcmpi(names, species), 1, 'first');
    if isempty(idx)
        name = '';
    else
        name = names{idx};
    end
end

function out = runCase(mesh, cfg, engine, chemistry, mode)
    area = sum(mesh.inletRingArea);
    rhoIn = engine.Pin / (cfg.Rgas * engine.Tin);
    u = max(engine.mdotAir / max(rhoIn * area, 1e-9), 0.5);

    dtAdv = cfg.cfl * mesh.dx / u;
    dtDiff = 0.24 * min(mesh.dx^2, mesh.dr^2) / ...
        max([cfg.DFuel, cfg.DO2, cfg.DNO, cfg.alpha]);
    dt = min([cfg.dtMax, dtAdv, dtDiff]);
    nSteps = max(1, ceil(cfg.tEnd / dt));
    dt = cfg.tEnd / nSteps;

    AFRst = 14.94; % c12h26 + 17.75 O2
    mdotFuelTarget = engine.mdotAir * cfg.phi / AFRst;

    Yfuel = zeros(mesh.nr, mesh.nx);
    YO2 = cfg.YO2In * ones(mesh.nr, mesh.nx);
    T = engine.Tin * ones(mesh.nr, mesh.nx);
    P = engine.Pin * ones(mesh.nr, mesh.nx);
    YNO = zeros(mesh.nr, mesh.nx);

    fuelConsumedField = zeros(mesh.nr, mesh.nx);

    singleInletFuel = zeros(mesh.nr, 1);
    sourceFuel = zeros(mesh.nr, mesh.nx);

    if strcmp(mode, 'single')
        singleInletFuel = buildSingleInletProfile(mesh, cfg, rhoIn, u, mdotFuelTarget);
        mdotFuelIn = rhoIn * u * sum(mesh.inletRingArea .* singleInletFuel);
    elseif strcmp(mode, 'dual')
        sourceFuel = buildDualSourceMap(mesh, cfg, rhoIn, mdotFuelTarget);
        mdotFuelIn = mdotFuelTarget;
    else
        error('Unknown mode: %s', mode);
    end

    for n = 1:nSteps
        rho = max(P ./ (cfg.Rgas * max(T, cfg.minTemperature)), 1e-6);

        % Inlet boundary conditions (x = 1)
        T(:, 1) = engine.Tin;
        P(:, 1) = engine.Pin;
        YO2(:, 1) = cfg.YO2In;
        if strcmp(mode, 'single')
            Yfuel(:, 1) = singleInletFuel;
        else
            Yfuel(:, 1) = 0.0;
        end

        [omegaFuel, omegaNO] = chemistryRates(T, Yfuel, YO2, P, chemistry, cfg);

        dYfdx = upstreamGradient(Yfuel, mesh.dx);
        dYO2dx = upstreamGradient(YO2, mesh.dx);
        dTdx = upstreamGradient(T, mesh.dx);
        dYNOdx = upstreamGradient(YNO, mesh.dx);

        lapYf = laplacianField(Yfuel, mesh.dr, mesh.dx);
        lapYO2 = laplacianField(YO2, mesh.dr, mesh.dx);
        lapT = laplacianField(T, mesh.dr, mesh.dx);
        lapYNO = laplacianField(YNO, mesh.dr, mesh.dx);

        if strcmp(mode, 'single')
            injectionSource = zeros(mesh.nr, mesh.nx);
        else
            injectionSource = sourceFuel ./ rho;
        end

        Yfuel = Yfuel + dt * (-u * dYfdx + cfg.DFuel * lapYf + injectionSource - omegaFuel);
        YO2 = YO2 + dt * (-u * dYO2dx + cfg.DO2 * lapYO2 - chemistry.fallbackStoichO2 * omegaFuel);
        T = T + dt * (-u * dTdx + cfg.alpha * lapT + (cfg.LHV / cfg.cp) * omegaFuel);
        YNO = YNO + dt * (-u * dYNOdx + cfg.DNO * lapYNO + omegaNO);

        % Bounds and positivity safeguards
        Yfuel = min(max(Yfuel, 0.0), cfg.maxMassFraction);
        YO2 = min(max(YO2, 0.0), cfg.maxMassFraction);
        YNO = min(max(YNO, 0.0), cfg.maxMassFraction);
        T = min(max(T, cfg.minTemperature), cfg.maxTemperature);

        fuelConsumedField = fuelConsumedField + dt * max(omegaFuel, 0.0);

        % March pressure field from inlet to outlet with thermal contribution
        for i = 2:mesh.nx
            deltaT = max(T(:, i-1) - engine.Tin, 0.0);
            localDrop = (cfg.pressureLossK * rho(:, i-1) * u^2 + ...
                cfg.pressureHeatK * engine.Pin * deltaT / (engine.Tin * cfg.length)) * mesh.dx;
            P(:, i) = max(cfg.pressureFloor, P(:, i-1) - localDrop);
        end
    end

    rhoOut = max(P(:, end) ./ (cfg.Rgas * max(T(:, end), cfg.minTemperature)), 1e-6);
    mdotOutletFuel = u * sum(mesh.inletRingArea .* rhoOut .* Yfuel(:, end));
    mdotFuelBurned = max(mdotFuelIn - mdotOutletFuel, 0.0);

    out = struct();
    out.case = mode;
    out.meshShape = [mesh.nr, mesh.nx];
    out.timeStep = dt;
    out.nSteps = nSteps;
    out.mdotFuelIn = mdotFuelIn;
    out.mdotFuelBurned = mdotFuelBurned;
    out.fuelConsumed = mdotFuelBurned * cfg.tEnd;
    out.outletTemperatureMean = mean(T(:, end));
    out.outletPressureMean = mean(P(:, end));
    out.outletTemperatureProfile = T(:, end);
    out.outletPressureProfile = P(:, end);
    out.outletFuelProfile = Yfuel(:, end);
    out.Yfuel = Yfuel;
    out.T = T;
    out.P = P;
    out.YNO = YNO;
    out.fuelConsumedField = fuelConsumedField;
    out.chemistryMode = chemistry.mode;
    out.hasMechanismNO = chemistry.hasMechanismNO;
end

function y = buildSingleInletProfile(mesh, cfg, rhoIn, u, mdotFuel)
    y = zeros(mesh.nr, 1);

    rFrac = ((1:mesh.nr)' - 0.5) / mesh.nr;
    mask = (rFrac >= cfg.singleRegion.rMinFrac) & (rFrac <= cfg.singleRegion.rMaxFrac);

    if ~any(mask)
        error('Single inlet region mask is empty.');
    end

    ringA = mesh.inletRingArea(mask);
    w = ringA / sum(ringA);
    y(mask) = mdotFuel * w ./ max(rhoIn * u * ringA, 1e-12);
    y = min(max(y, 0.0), 0.35);
end

function source = buildDualSourceMap(mesh, cfg, rhoIn, mdotFuel)
    source = zeros(mesh.nr, mesh.nx);

    for k = 1:2
        region = cfg.dualRegions{k};
        xFrac = (mesh.X / (mesh.nx * mesh.dx));
        rFrac = (mesh.R / (mesh.nr * mesh.dr));

        mask = (xFrac >= region.xMinFrac) & (xFrac <= region.xMaxFrac) & ...
               (rFrac >= region.rMinFrac) & (rFrac <= region.rMaxFrac);

        if ~any(mask(:))
            error('Dual inlet region %d mask is empty.', k);
        end

        vol = mesh.cellVolume(mask);
        mdotRegion = cfg.dualSplit(k) * mdotFuel;
        source(mask) = source(mask) + mdotRegion / max(sum(vol), 1e-12);
    end

    % Convert from volumetric source to conservative bounded value
    source = min(max(source, 0.0), 20.0 * rhoIn);
end

function grad = upstreamGradient(field, dx)
    grad = zeros(size(field));
    grad(:, 2:end) = (field(:, 2:end) - field(:, 1:end-1)) / dx;
    grad(:, 1) = 0.0;
end

function lap = laplacianField(field, dr, dx)
    frp = [field(2:end, :); field(end, :)];
    frm = [field(1, :); field(1:end-1, :)];
    fzp = [field(:, 2:end), field(:, end)];
    fzm = [field(:, 1), field(:, 1:end-1)];

    lap = (frp - 2.0 * field + frm) / dr^2 + (fzp - 2.0 * field + fzm) / dx^2;
end

function [omegaFuel, omegaNO] = chemistryRates(T, Yfuel, YO2, P, chemistry, cfg)
    omegaFuel = zeros(size(T));
    omegaNO = zeros(size(T));

    if chemistry.useCantera
        [omegaFuel, omegaNO] = canteraLookupRates(T, Yfuel, YO2, P, chemistry, cfg);
        return;
    end

    % Fallback bounded Arrhenius source
    arr = chemistry.fallbackA .* exp(-chemistry.fallbackE ./ (8.314 * max(T, 300.0)));
    omegaFuel = arr .* Yfuel .* min(YO2 / cfg.YO2In, 1.0);
    omegaFuel = min(max(omegaFuel, 0.0), cfg.maxFuelRate);

    if any(T(:) > 1800.0)
        omegaNO = thermalNORate(T, YO2);
    end
end

function [omegaFuel, omegaNO] = canteraLookupRates(T, Yfuel, YO2, P, chemistry, cfg)
    omegaFuel = zeros(size(T));
    omegaNO = zeros(size(T));

    fuelName = chemistry.fuelName;
    o2Name = chemistry.o2Name;
    n2Name = chemistry.n2Name;

    for i = 1:numel(T)
        Ti = max(T(i), 300.0);
        Pi = max(P(i), cfg.pressureFloor);
        Yfi = max(Yfuel(i), 0.0);
        YO2i = max(YO2(i), 0.0);
        YN2i = max(1.0 - Yfi - YO2i, 1e-12);

        key = sprintf('%s|%s|%d|%d|%d|%d', chemistry.mechanismPath, chemistry.phase, ...
            round(Ti), round(1e4 * Yfi), round(1e4 * YO2i), round(Pi / 100));

        cached = canteraLookupCache(key);
        if ~isempty(cached)
            omegaFuel(i) = cached(1);
            omegaNO(i) = cached(2);
            continue;
        end

        xFuel = Yfi / max(chemistry.MWFuel, 1e-12);
        xO2 = YO2i / max(chemistry.MWO2, 1e-12);
        xN2 = YN2i / max(chemistry.MWN2, 1e-12);
        xSum = xFuel + xO2 + xN2;

        comp = sprintf('%s:%0.12g, %s:%0.12g, %s:%0.12g', ...
            fuelName, xFuel / xSum, o2Name, xO2 / xSum, n2Name, xN2 / xSum);

        try
            chemistry.gas.TPX = {Ti, Pi, comp};
            rho = max(chemistry.gas.density, 1e-6);
            wdot = chemistry.gas.netProdRates;

            wFuel = max(-wdot(chemistry.fuelIndex) * chemistry.MWFuel / rho, 0.0);
            wFuel = min(wFuel, cfg.maxFuelRate);

            if chemistry.hasMechanismNO
                wNO = wdot(chemistry.noIndex) * chemistry.MWNO / rho;
            else
                wNO = thermalNORate(Ti, YO2i);
            end

            omegaFuel(i) = wFuel;
            omegaNO(i) = wNO;
            canteraLookupCache(key, [wFuel, wNO]);

        catch
            arr = chemistry.fallbackA * exp(-chemistry.fallbackE / (8.314 * Ti));
            omegaFuel(i) = min(arr * Yfi * min(YO2i / cfg.YO2In, 1.0), cfg.maxFuelRate);
            omegaNO(i) = thermalNORate(Ti, YO2i);
        end
    end
end

function value = canteraLookupCache(key, varargin)
    persistent CACHE
    if isempty(CACHE)
        CACHE = containers.Map('KeyType', 'char', 'ValueType', 'any');
    end

    if nargin == 2
        CACHE(key) = varargin{1};
        value = [];
        return;
    end

    if isKey(CACHE, key)
        value = CACHE(key);
    else
        value = [];
    end
end

function rate = thermalNORate(T, YO2)
    % Thermal-NO estimate (fallback) for trend visualization only.
    R = 8.314;
    Tloc = max(T, 300.0);
    rate = 2.0e-3 .* YO2 .* exp(-7.0e4 ./ (R .* Tloc));
end

function plotCaseFields(results, mesh, engine, cfg)
    single = results.single;
    dual = results.dual;

    rNorm = ((1:mesh.nr) - 0.5) / mesh.nr;
    xNorm = ((1:mesh.nx) - 0.5) / mesh.nx;

    figure('Color', 'w', 'Position', [100 100 1400 850]);

    subplot(2, 3, 1)
    imagesc(xNorm, rNorm, single.fuelConsumedField);
    axis xy;
    title('Single inlet: fuel consumption field');
    xlabel('Axial fraction');
    ylabel('Radial fraction');
    colorbar;

    subplot(2, 3, 2)
    imagesc(xNorm, rNorm, single.T);
    axis xy;
    title('Single inlet: temperature field [K]');
    xlabel('Axial fraction');
    ylabel('Radial fraction');
    colorbar;

    subplot(2, 3, 3)
    imagesc(xNorm, rNorm, single.P / 1e5);
    axis xy;
    title('Single inlet: pressure field [bar]');
    xlabel('Axial fraction');
    ylabel('Radial fraction');
    colorbar;

    subplot(2, 3, 4)
    imagesc(xNorm, rNorm, dual.fuelConsumedField);
    axis xy;
    title('Dual inlet: fuel consumption field');
    xlabel('Axial fraction');
    ylabel('Radial fraction');
    colorbar;

    subplot(2, 3, 5)
    imagesc(xNorm, rNorm, dual.T);
    axis xy;
    title('Dual inlet: temperature field [K]');
    xlabel('Axial fraction');
    ylabel('Radial fraction');
    colorbar;

    subplot(2, 3, 6)
    imagesc(xNorm, rNorm, dual.P / 1e5);
    axis xy;
    title('Dual inlet: pressure field [bar]');
    xlabel('Axial fraction');
    ylabel('Radial fraction');
    colorbar;

    sgtitle(sprintf('%s | 2-D combustor model (%s chemistry)', ...
        engine.key, results.metadata.chemistryMode));

    figure('Color', 'w', 'Position', [120 120 1100 420]);

    subplot(1, 3, 1)
    bar([single.mdotFuelBurned, dual.mdotFuelBurned]);
    set(gca, 'XTickLabel', {'Single', 'Dual'});
    ylabel('kg/s');
    title('Fuel consumption rate');

    subplot(1, 3, 2)
    plot(rNorm, single.outletTemperatureProfile, 'b-', 'LineWidth', 1.8);
    hold on;
    plot(rNorm, dual.outletTemperatureProfile, 'r--', 'LineWidth', 1.8);
    grid on;
    xlabel('Radial fraction');
    ylabel('Outlet temperature [K]');
    title('Outlet temperature profile');
    legend('Single', 'Dual', 'Location', 'best');

    subplot(1, 3, 3)
    plot(rNorm, single.outletPressureProfile / 1e5, 'b-', 'LineWidth', 1.8);
    hold on;
    plot(rNorm, dual.outletPressureProfile / 1e5, 'r--', 'LineWidth', 1.8);
    grid on;
    xlabel('Radial fraction');
    ylabel('Outlet pressure [bar]');
    title('Outlet pressure profile');
    legend('Single', 'Dual', 'Location', 'best');

    sgtitle(sprintf('Assumed engine inlet: T = %.0f K, P = %.1f bar | \phi = %.2f', ...
        engine.Tin, engine.Pin / 1e5, cfg.phi));
end
