%% SETUP  Add ShepherdLab LW FEA-OPT source directories to the MATLAB path.
%
%   Run this script once per session (or add it to your startup.m).
%
%       >> run('path/to/lw-feaopt/setup.m')
%
%   Adds:
%     src/        — core classes and optimiser
%     io/         — geometry importers and profile exporter
%     examples/   — demo scripts
%
%   See also: demo_beam_optimization

root = fileparts(mfilename('fullpath'));

addpath(fullfile(root, 'src'));
addpath(fullfile(root, 'io'));
addpath(fullfile(root, 'examples'));

% Ensure profiles output directory exists
profiles_dir = fullfile(root, 'profiles');
if ~exist(profiles_dir, 'dir')
    mkdir(profiles_dir);
end

fprintf('ShepherdLab LW FEA-OPT paths added.\n');
fprintf('  src/        — core library\n');
fprintf('  io/         — importers & profile exporter\n');
fprintf('  examples/   — demo scripts\n');
fprintf('  profiles/   — exported CAD profiles\n');
fprintf('\nRun "demo_beam_optimization" to get started.\n');