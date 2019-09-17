%% Setup VLFeat toolbox and add other papers' codes.
addpath('modelspecific'); addpath('multigs');  % for feature match and homography

addpath('texture_mapping'); % for our texture mapping

addpath('LSD_matlab'); addpath('LineMatching'); % for line segments detection and match

% Setup VLFeat toolbox
run('vlfeat-0.9.21/toolbox/vl_setup');