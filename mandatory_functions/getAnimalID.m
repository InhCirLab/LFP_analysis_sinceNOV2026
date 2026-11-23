function id = getAnimalID(name)
% getAnimalID  Extract the animal ID from a session folder name.
%
%   id = getAnimalID(name) returns the last dash-delimited chunk of a
%   session name. This matches the lab's naming convention where the animal
%   ID is always the final segment after the last '-'.
%
%   Examples:
%       getAnimalID('2026-06-09_16-47-11-F016')  ->  'F016'
%       getAnimalID('2025-11-12-124020-E165')    ->  'E165'
%
%   Input:
%       name  char or string, a single session/folder name
%   Output:
%       id    char, the trailing segment after the final '-'
%
%   Note: returns the whole name (trimmed) if there is no '-' in it.

    name  = char(name);
    parts = strsplit(name, '-');
    id    = strtrim(parts{end});
end
