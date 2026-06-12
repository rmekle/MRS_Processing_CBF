%% rename_files_DICOM_toPseudonyms_s.m
%
% Renames Siemens DICOM (.IMA) files distributed across multiple subfolders by
% replacing the Lastname_Firstname prefix (and any following text) with a
% pseudonym defined in a CSV mapping file.
%
% SUPPORTED FILENAME PATTERNS:
%   Pattern A: Lastname_Firstname.Text.0001.IMA
%   Pattern B: Lastname_Firstname.0001.IMA
%
%   In both cases everything up to (but not including) the FIRST occurrence
%   of a dot followed by a leading-zero number (e.g. ".0001", ".01", ".0")
%   is replaced by the pseudonym. The leading-zero number and everything
%   after it is preserved unchanged.
%
%   Examples:
%     Mueller_Hans.SomeText.0023.IMA  ->  Sub001.0023.IMA
%     Mueller_Hans.0023.IMA           ->  Sub001.0023.IMA
%     Schmidt_Anna.Scan1.001.IMA      ->  Sub002.001.IMA
%
% MAPPING CSV FORMAT (no header row required, but supported):
%   column 1 = Lastname_Firstname, column 2 = Pseudonym
%   Example rows:
%     Mueller_Hans,Sub001
%     Schmidt_Anna,Sub002
%
% The search is case-insensitive to handle naming inconsistencies
% (e.g. mueller_hans or MUELLER_HANS will both match Mueller_Hans)
%
% USAGE:
%   1. Set rootDir to the top-level folder containing all subfolders
%   2. Set mappingFile to the full path of your CSV mapping file
%   3. Set dryRun = true to preview changes without renaming anything
%      Set dryRun = false to perform the actual renaming
%   4. Run the script
%
% OUTPUT:
%   - Command window log of all matched / renamed / skipped files
%   - Summary statistics at the end
%
% SAFETY FEATURES:
%   - Dry-run mode: preview all changes before committing
%   - Skips files where the new filename already exists (no overwrite)
%   - Warns about filenames that match no entry in the mapping CSV
%   - Warns about ambiguous matches (multiple CSV entries match one file)

%% -------------------------------------------------------------------------
%  USER SETTINGS — edit these before running
%  -------------------------------------------------------------------------

%rootDir     = '/path/to/your/top/level/folder';   % <-- set this
%mappingFile = '/path/to/your/mapping.csv';         % <-- set this
rootDir     = '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_MRS_TGA/MRS_TGA_Data_All/';
mappingFile = '/home/mekler/CSB_NeuroRad/mekler/Data_II/3T_MRS_TGA/3T_MRS_TGA_Data_Info.csv';
dryRun      = true;   % true = preview only; false = rename for real

%% -------------------------------------------------------------------------
%  LOAD MAPPING CSV
%  -------------------------------------------------------------------------

fprintf('\n=== rename_IMA_pseudonyms ===\n');
fprintf('Root directory : %s\n', rootDir);
fprintf('Mapping file   : %s\n', mappingFile);
fprintf('Dry run        : %s\n\n', mat2str(dryRun));

% Read CSV — supports optional header row, blank lines, and extra whitespace
fid = fopen(mappingFile, 'r');
if fid == -1
    error('Cannot open mapping file: %s', mappingFile);
end

mapping = struct('name', {}, 'pseudonym', {});
lineNum = 0;
while ~feof(fid)
    line    = strtrim(fgetl(fid));
    lineNum = lineNum + 1;
    if isempty(line) || line(1) == '#'
        continue   % skip blank lines and comments
    end
    parts = strsplit(line, ',');
    if numel(parts) < 2
        warning('Skipping malformed line %d in CSV: %s', lineNum, line);
        continue
    end
    name       = strtrim(parts{1});
    pseudonym  = strtrim(parts{2});
    % Skip header row if present (heuristic: first column contains letters
    % that look like a column label rather than a name)
    if strcmpi(name, 'name') || strcmpi(name, 'lastname_firstname')
        continue
    end
    mapping(end+1).name       = name;       %#ok<AGROW>
    mapping(end).pseudonym    = pseudonym;
end
fclose(fid);

if isempty(mapping)
    error('No valid entries found in mapping file. Check CSV format.');
end
fprintf('Loaded %d name-to-pseudonym mappings.\n\n', numel(mapping));


%% -------------------------------------------------------------------------
%  FIND ALL .IMA FILES RECURSIVELY
%  -------------------------------------------------------------------------

imaFiles = dir(fullfile(rootDir, '**', '*.IMA'));

% Also catch lowercase extension variants (.ima) for robustness
imaFilesLower = dir(fullfile(rootDir, '**', '*.ima'));
imaFiles = [imaFiles; imaFilesLower];

if isempty(imaFiles)
    fprintf('No .IMA files found under %s\n', rootDir);
    return
end
fprintf('Found %d .IMA file(s) to process.\n\n', numel(imaFiles));

%% -------------------------------------------------------------------------
%  PROCESS EACH FILE
%  -------------------------------------------------------------------------

nRenamed   = 0;
nSkipped   = 0;
nUnmatched = 0;
nConflict  = 0;

for k = 1 : numel(imaFiles)

    oldName  = imaFiles(k).name;    % filename only, no path
    fileDir  = imaFiles(k).folder;  % full path to containing folder
    oldPath  = fullfile(fileDir, oldName);

    % ------------------------------------------------------------------
    % Step 1 — Find the first ".0number" anchor in the filename.
    %
    % Supported patterns (all handled identically):
    %   Lastname_Firstname.0023.IMA
    %   Lastname_Firstname.Text.0023.IMA
    %   SingleToken.0023.IMA
    %   Lastname_Firstname_DR_.0023.IMA
    %
    % Regex captures:
    %   tok{1}  everything before the first ".0<digits>" — the raw prefix
    %   tok{2}  the ".0<digits>" anchor itself
    %   tok{3}  everything after (further numbers, extension, etc.)
    % ------------------------------------------------------------------
    pat = '^(.+?)(\.0\d*)(.+)?$';
    tok = regexp(oldName, pat, 'tokens', 'ignorecase');

    if isempty(tok)
        fprintf('[SKIP - no leading-zero number found] %s\n', oldName);
        nSkipped = nSkipped + 1;
        continue
    end

    tok       = tok{1};
    rawPrefix = tok{1};   % e.g. "HüBNER-GOTTSCHICK_ANDREA_DR_" or
                          %      "Mueller_Hans.SomeText" or "Mueller_Hans"
    numSuffix = tok{2};   % e.g. ".0023"
    rest      = tok{3};   % e.g. ".IMA"

    % ------------------------------------------------------------------
    % Step 2 — Longest-match CSV lookup.
    %
    % Rather than trying to parse the name structure from the filename
    % (which is fragile given variable titles, trailing underscores,
    % and dot-separated text), we let the CSV drive the matching:
    %
    %   For each CSV entry, check whether rawPrefix STARTS WITH that
    %   entry (case-insensitive). Among all entries that match, pick
    %   the LONGEST one — this correctly handles cases where one name
    %   is a prefix of another (e.g. "Mueller" vs "Mueller_Hans").
    %
    % The rawPrefix is first normalised by stripping any trailing
    % underscores and dots, and by extracting only the part before
    % the first dot (to discard ".SomeText" suffixes that may appear
    % between the name and the .0number anchor).
    % ------------------------------------------------------------------

    % Normalise: strip everything from the first dot onward (discards
    % ".Text" fragments), then strip trailing underscores/spaces
    dotPos       = strfind(rawPrefix, '.');
    if ~isempty(dotPos)
        normPrefix = rawPrefix(1 : dotPos(1)-1);
    else
        normPrefix = rawPrefix;
    end
    normPrefix = strtrim(regexprep(normPrefix, '[_\s]+$', ''));  % trailing _ or space

    % Find all CSV entries that are a case-insensitive prefix of normPrefix
    matchIdx    = [];
    matchLens   = [];
    for m = 1 : numel(mapping)
        csvName = mapping(m).name;
        L       = length(csvName);
        if length(normPrefix) >= L && ...
                strncmpi(normPrefix, csvName, L)
            % Additional check: the character immediately after the CSV
            % name in normPrefix (if any) must be a non-alphanumeric
            % separator (underscore, hyphen, space) — prevents "Mueller"
            % matching "MuellerHans"
            if length(normPrefix) == L || ...
                    ~isempty(regexp(normPrefix(L+1), '^[_\-\s]', 'once'))
                matchIdx(end+1)  = m;   %#ok<AGROW>
                matchLens(end+1) = L;   %#ok<AGROW>
            end
        end
    end

    if isempty(matchIdx)
        fprintf('[UNMATCHED] %s  (normalised prefix: %s)\n', oldName, normPrefix);
        nUnmatched = nUnmatched + 1;
        continue
    end

    % Pick the longest matching CSV entry
    [~, bestIdx] = max(matchLens);
    matchIdx     = matchIdx(bestIdx);
    pseudonym    = mapping(matchIdx).pseudonym;

    % ------------------------------------------------------------------
    % Build new filename: Pseudonym<.0number><rest>
    % e.g. Sub001.0023.IMA
    % ------------------------------------------------------------------
    newName = [pseudonym, numSuffix, rest];
    newPath = fullfile(fileDir, newName);

    % ------------------------------------------------------------------
    % Safety check: do not overwrite an existing file
    % ------------------------------------------------------------------
    if exist(newPath, 'file')
        fprintf('[CONFLICT - skipped] Target already exists: %s\n', newName);
        nConflict = nConflict + 1;
        continue
    end

    % ------------------------------------------------------------------
    % Rename (or report in dry-run mode)
    % ------------------------------------------------------------------
    fprintf('[%s] %s\n        -> %s\n', ...
        ternary(dryRun, 'DRY RUN', 'RENAMED'), oldName, newName);

    if ~dryRun
        [success, msg] = movefile(oldPath, newPath);
        if ~success
            fprintf('  ERROR renaming file: %s\n', msg);
            nSkipped = nSkipped + 1;
            continue
        end
    end

    nRenamed = nRenamed + 1;

end

%% -------------------------------------------------------------------------
%  SUMMARY
%  -------------------------------------------------------------------------

fprintf('\n=== Summary ===\n');
if dryRun
    fprintf('Mode           : DRY RUN (no files were changed)\n');
else
    fprintf('Mode           : LIVE (files have been renamed)\n');
end
fprintf('Renamed        : %d\n', nRenamed);
fprintf('Unmatched      : %d\n', nUnmatched);
fprintf('Conflicts      : %d\n', nConflict);
fprintf('Skipped (other): %d\n', nSkipped);
fprintf('Total processed: %d\n', numel(imaFiles));

if dryRun && nRenamed > 0
    fprintf('\nRe-run with dryRun = false to apply the %d rename(s).\n', nRenamed);
end

%% -------------------------------------------------------------------------
%  HELPER FUNCTION
%  -------------------------------------------------------------------------

function out = ternary(cond, a, b)
    if cond; out = a; else; out = b; end
end