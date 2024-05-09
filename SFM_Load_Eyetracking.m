% Load in eyetracking data for SFM
% KWK - 20201030

function [ETData] = SFM_Load_Eyetracking


clear all; close all;

% Make the subject list out of the file names
options.curDur = pwd;
% Find and load in the eye tracking data
cd /home/shaw-raid1/data/psychophysics/eyetracking/
fileNames = dir('*');

cd(options.curDur)








end


