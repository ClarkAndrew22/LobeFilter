% Andrew Clark
% Discription: Code to process fluorescent particle through Manta Ray Filter and store
% particle counts in csv for statistical analysis in another code and/or JMP

% Experiment: Device 1 with 25 and 15 um particle collection 
% Experiment Date: 07-05-2021
% Written: 07-05-2021


% Transfer between folders
run = {'R1','R2','R3'};

for kk = 1:3
    % Change directory to experiment
    cd(run{kk})
    
    close all;
%     clear;
%     clc;

    % Start of new written code %
    Inlet_Matrix_red = zeros(3,1);
    Inlet_Matrix_green = zeros(3,1);
    out1_matrix_red = zeros(3,1);
    out1_matrix_green = zeros(3,1);
    out2_matrix_red = zeros(3,1);
    out2_matrix_green = zeros(3,1);

    for jj = 1:3
        inlet_im = imread(['8-in-', num2str(jj) '.tif']);
        inlet_red = inlet_im(:,:,1);
        inlet_red_bw = imbinarize(inlet_red,.3);
        inlet_red_bw = bwareaopen(inlet_red_bw,3);
        inlet_green = inlet_im(:,:,2);
        inlet_green_bw = imbinarize(inlet_green,.3);
        inlet_green_bw = bwareaopen(inlet_green_bw,3);

        inlet_red_count = length(regionprops(inlet_red_bw,'Area'));
        inlet_green_count = length(regionprops(inlet_green_bw,'Area'));

        Inlet_Matrix_red(jj) = inlet_red_count;
        Inlet_Matrix_green(jj) = inlet_green_count;

        out1_im = imread(['8-1-' num2str(jj) '.tif']); 
        out2_im = imread(['8-2-' num2str(jj) '.tif']); 
        out1_red = out1_im(:,:,1);
        out2_red = out2_im(:,:,1);
        out1_red_bw = imbinarize(out1_red,.3);
        out1_red_bw = bwareaopen(out1_red_bw,3);
        out2_red_bw = imbinarize(out2_red,.3);
        out2_red_bw = bwareaopen(out2_red_bw,3);
        out1_green = out1_im(:,:,2);
        out2_green = out2_im(:,:,2);
        out1_green_bw = imbinarize(out1_green,.3);
        out1_green_bw = bwareaopen(out1_green_bw,3);
        out2_green_bw = imbinarize(out2_green,.3);
        out2_green_bw = bwareaopen(out2_green_bw,3);

        out1_red_count = length(regionprops(out1_red_bw));
        out2_red_count = length(regionprops(out2_red_bw));
        out1_green_count = length(regionprops(out1_green_bw));
        out2_green_count = length(regionprops(out2_green_bw));

        out1_matrix_red(jj) = out1_red_count;
        out1_matrix_green(jj) = out1_green_count;
        out2_matrix_red(jj) = out2_red_count;
        out2_matrix_green(jj) = out2_green_count;
    end

    % Collect
    inRed(kk,1) = mean(Inlet_Matrix_red);
    inGreen(kk,1) = mean(Inlet_Matrix_green);
    o1red(kk,1) = mean(out1_matrix_red);
    o1green(kk,1) = mean(out1_matrix_green);
    o2red(kk,1) = mean(out2_matrix_red);
    o2green(kk,1) = mean(out2_matrix_green);
    
    inRedstd(kk,1) = std(Inlet_Matrix_red);
    inGreenstd(kk,1) = std(Inlet_Matrix_green);
    o1redstd(kk,1) = std(out1_matrix_red);
    o1greenstd(kk,1) = std(out1_matrix_green);
    o2redstd(kk,1) = std(out2_matrix_red);
    o2greenstd(kk,1) = std(out2_matrix_green);

    % Get efficiency and concentration ratio
    redEff(kk) = (1-o2red(kk)./inRed(kk)).*100;
    greenEff(kk) = (1-o2green(kk)./inGreen(kk)).*100;
    redCR(kk) = o1red(kk)./inRed(kk);
    greenCR(kk) = o1green(kk)./inGreen(kk);
    rPur(kk) = o1red(kk)./(o1red(kk)+o1green(kk));
    gPur(kk) = o2green(kk)./(o2green(kk)+o2red(kk));
%     rSE(kk) = (o1red./o1green)./(inRed./inGreen);
%     gSE(kk) = (o2green./o2red)./(inGreen./inRed);
    % Change back to main folder
    cd ..

end

%% Gather data
redEffMean = mean(redEff);
greenEffMean = mean(greenEff);
redEffStd = std(redEff);
greenEffStd = std(greenEff);

rCRmean = mean(redCR);
gCRmean = mean(greenCR);
rCRstd = std(redCR);
gCRstd = std(greenCR);

rPmean = mean(rPur);
gPmean = mean(gPur);
rPstd = std(rPur);
gPstd = std(gPur);

% redSE = mean(rSE);
% redSEstd = std(rSE);
% greenSE = mean(gSE);
% greenSEstd = std(gSE);

% Store
% save('runData','redEffMean','greenEffMean','redEffStd','greenEffStd','rCRmean','gCRmean','rCRstd','gCRstd');





















