%% Community metabolic model reconstruction and flux export
% This script documents how sample-specific gut microbiome community models
% were reconstructed from detected species and normalized abundance data.
%
% Requirements:
%   - MATLAB
%   - COBRA Toolbox
%   - mgPipe / Microbiome Modeling Toolbox
%   - APOLLO gut microbiome genome-scale reconstruction catalogue
%
% Edit the configuration block below before running.

%% ==============================
% 0. Configuration
%% ==============================

speciesListFile = fullfile('data', 'processed', 'species_list.tsv');
apolloDir       = fullfile('external', 'APOLLO_models');
modelOutDir     = fullfile('results', 'community_models', 'selected_species_models');
matchLogFile    = fullfile(modelOutDir, 'match_log.tsv');
mappingOutFile  = fullfile(modelOutDir, 'species_to_GEM_mapping.tsv');

abundanceFile   = fullfile('data', 'processed', 'NormCoverage_Cleaned.csv');
dietFile        = fullfile('external', 'diets', 'AverageEuropeanDiet.txt');
communityOutDir = fullfile('results', 'community_models', 'mgpipe_output');
simulationOutDir = fullfile('results', 'simulation_results');

numWorkers = 4;

if ~exist(modelOutDir, 'dir'); mkdir(modelOutDir); end
if ~exist(communityOutDir, 'dir'); mkdir(communityOutDir); end
if ~exist(simulationOutDir, 'dir'); mkdir(simulationOutDir); end

%% ==============================
% 1. Match detected species to APOLLO models
%% ==============================
% species_list.tsv is expected to contain either:
%   Genus<TAB>Species
% or a single species binomial per row.

fid = fopen(speciesListFile);
if fid < 0
    error('Could not open species list: %s', speciesListFile);
end
raw = textscan(fid, '%s %s', 'Delimiter', '\t');
fclose(fid);

if isempty(raw{2})
    speciesList = strtrim(raw{1});
else
    speciesList = strcat(strtrim(raw{1}), {' '}, strtrim(raw{2}));
end
speciesList = speciesList(~cellfun('isempty', speciesList));
speciesList = unique(speciesList);

fid = fopen(matchLogFile, 'w');
fprintf(fid, 'MAG_Species\tModel_Selected\tMatch_Level\tNotes\n');

for i = 1:length(speciesList)
    taxon = strtrim(speciesList{i});
    fprintf('Processing species: %s\n', taxon);

    tokenized = split(taxon);
    if numel(tokenized) < 2
        fprintf(fid, '%s\tNA\tunmatched\tinvalid species name\n', taxon);
        continue
    end

    genus = tokenized{1};
    species = tokenized{2};

    % First try exact species-level APOLLO pan-model match.
    speciesPattern = sprintf('pan%s_%s*.mat', genus, species);
    speciesMatch = dir(fullfile(apolloDir, speciesPattern));

    if ~isempty(speciesMatch)
        copyfile(fullfile(apolloDir, speciesMatch(1).name), modelOutDir);
        fprintf(fid, '%s\t%s\tspecies\texact match\n', taxon, speciesMatch(1).name);
        continue
    end

    % If no species-level model exists, use the first available genus-level pan-model.
    genusPattern = sprintf('pan%s_*.mat', genus);
    genusMatch = dir(fullfile(apolloDir, genusPattern));

    if ~isempty(genusMatch)
        copyfile(fullfile(apolloDir, genusMatch(1).name), modelOutDir);
        fprintf(fid, '%s\t%s\tgenus\tfallback genus match\n', taxon, genusMatch(1).name);
        continue
    end

    % Optional fallback: attempt pan-genus model construction.
    fprintf('No existing model found for %s. Attempting pan-genus build.\n', taxon);
    try
        try
            createPanModels(apolloDir, modelOutDir, 'genus', 'local', numWorkers);
        catch
            createPanModels(apolloDir, modelOutDir, 'genus', 'local');
        end
        fprintf(fid, '%s\t%s\tpan-genus-created\tnew build attempted\n', taxon, genus);
    catch ME
        warning('Pan-genus build failed for %s: %s', genus, ME.message);
        fprintf(fid, '%s\tNA\tunmatched\tcreatePanModels failed: %s\n', taxon, ME.message);
    end
end

fclose(fid);
fprintf('Model matching log written to: %s\n', matchLogFile);

%% ==============================
% 2. Save species-to-GEM mapping table
%% ==============================

matchTbl = readtable(matchLogFile, 'FileType', 'text', 'Delimiter', '\t');
requiredCols = {'MAG_Species', 'Model_Selected'};
if ~all(ismember(requiredCols, matchTbl.Properties.VariableNames))
    error('match_log.tsv must contain MAG_Species and Model_Selected columns.');
end

mapTbl = matchTbl(:, requiredCols);
[~, ia] = unique(mapTbl.MAG_Species, 'stable');
mapTbl = mapTbl(ia, :);
writetable(mapTbl, mappingOutFile, 'FileType', 'text', 'Delimiter', '\t');
fprintf('Species-to-GEM mapping saved to: %s\n', mappingOutFile);

%% ==============================
% 3. Reconstruct sample-specific communities with mgPipe
%% ==============================
% abundanceFile should contain normalized species abundance/coverage values
% with model-compatible species identifiers.

initMgPipe(modelOutDir, abundanceFile, true, ...
    'NumWorkers', numWorkers, ...
    'rDiet', true, ...
    'resPath', communityOutDir, ...
    'DietFilePath', dietFile);

fprintf('mgPipe community reconstruction complete.\n');

%% ==============================
% 4. Run FBA on reconstructed community models and export flux summaries
%% ==============================

modelFiles = dir(fullfile(communityOutDir, '*.mat'));
if isempty(modelFiles)
    warning('No community model .mat files found in %s', communityOutDir);
end

dietData = readtable(dietFile, 'Delimiter', '\t', 'ReadVariableNames', false);
dietMetabolites = dietData.Var1;
dietAmounts = dietData.Var2;

allMetaboliteFluxes = table();
allReactionFluxes = table();

for k = 1:length(modelFiles)
    modelFileName = modelFiles(k).name;
    modelFilePath = fullfile(communityOutDir, modelFileName);
    fprintf('Loading model: %s\n', modelFileName);

    loaded = load(modelFilePath);
    if isfield(loaded, 'model')
        model = loaded.model;
    else
        warning('No variable named model found in %s. Skipping.', modelFileName);
        continue
    end

    % Apply diet constraints to matching exchange reactions.
    for i = 1:length(dietMetabolites)
        metID = dietMetabolites{i};
        amount = dietAmounts(i);
        reactionIdx = find(contains(model.rxns, metID));
        if ~isempty(reactionIdx)
            model.lb(reactionIdx) = -amount;
            model.ub(reactionIdx) = 0;
        end
    end

    fprintf('Running FBA optimization for %s\n', modelFileName);
    solution = optimizeCbModel(model);

    if isempty(solution.x)
        warning('Optimization failed for model: %s', modelFileName);
        continue
    end

    metaboliteFluxes = containers.Map();
    reactionFluxes = containers.Map();

    for i = 1:length(model.mets)
        metID = model.mets{i};
        [~, reactionIndices] = find(model.S(i, :));
        totalFlux = sum(solution.x(reactionIndices));
        metaboliteFluxes(metID) = totalFlux;
    end

    for i = 1:length(model.rxns)
        reactionFluxes(model.rxns{i}) = solution.x(i);
    end

    metaboliteIDs = keys(metaboliteFluxes);
    metaboliteValues = cell2mat(values(metaboliteFluxes));
    modelMetaboliteFluxTable = table(
        metaboliteIDs', metaboliteValues', ...
        'VariableNames', {'Metabolite', matlab.lang.makeValidName(modelFileName)});

    reactionIDs = keys(reactionFluxes);
    reactionValues = cell2mat(values(reactionFluxes));
    modelReactionFluxTable = table(
        reactionIDs', reactionValues', ...
        'VariableNames', {'Reaction', matlab.lang.makeValidName(modelFileName)});

    if isempty(allMetaboliteFluxes)
        allMetaboliteFluxes = modelMetaboliteFluxTable;
    else
        allMetaboliteFluxes = outerjoin(allMetaboliteFluxes, modelMetaboliteFluxTable, ...
            'Keys', 'Metabolite', 'MergeKeys', true);
    end

    if isempty(allReactionFluxes)
        allReactionFluxes = modelReactionFluxTable;
    else
        allReactionFluxes = outerjoin(allReactionFluxes, modelReactionFluxTable, ...
            'Keys', 'Reaction', 'MergeKeys', true);
    end
end

if ~isempty(allMetaboliteFluxes)
    numericColsMetabolites = varfun(@isnumeric, allMetaboliteFluxes, 'OutputFormat', 'uniform');
    allMetaboliteFluxes{:, numericColsMetabolites} = fillmissing(
        allMetaboliteFluxes{:, numericColsMetabolites}, 'constant', 0);
end

if ~isempty(allReactionFluxes)
    numericColsReactions = varfun(@isnumeric, allReactionFluxes, 'OutputFormat', 'uniform');
    allReactionFluxes{:, numericColsReactions} = fillmissing(
        allReactionFluxes{:, numericColsReactions}, 'constant', 0);
end

combinedMetaboliteSummaryFile = fullfile(simulationOutDir, 'Combined_Metabolite_Uptake_Production_Summary.csv');
combinedReactionSummaryFile = fullfile(simulationOutDir, 'Combined_Reaction_Flux_Summary.csv');

writetable(allMetaboliteFluxes, combinedMetaboliteSummaryFile);
writetable(allReactionFluxes, combinedReactionSummaryFile);

fprintf('Metabolite flux summary written to: %s\n', combinedMetaboliteSummaryFile);
fprintf('Reaction flux summary written to: %s\n', combinedReactionSummaryFile);
fprintf('Community reconstruction and simulation complete.\n');
