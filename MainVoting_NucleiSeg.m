%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is designed by Hongming Xu,
% Deptment of Eletrical and Computer Engineering,
% University of Alberta, Canada.  1th April, 2016
% If you have any problem feel free to contact me.
% Please address questions or comments to: mxu@ualberta.ca

% The techique is mainly based on the following paper:
% Xu, Hongming, et al. "An efficient technique for nuclei segmentation based on ellipse descriptor analysis and improved seed detection algorithm." (2014).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MainVoting_NucleiSeg
close all;
clear all;
%% %% add the function into MATLAB searching path and enter the test dataset
p = mfilename('fullpath');
t = findstr(p,'\');
p = p(1:t(end));
addpath(p);
cd(p);
cd('dataII');
list=dir('*.tif');


for i=1:length(list)
    filename=list(i).name;
    I=imread(filename);
    I2=normalizeStaining(I);
    R=I2(:,:,1);
    
    ac=30;
    figure,imshow(I);
    
    %% Image binarization -- segment nuclei regions into foreground
    bwf=~im2bw(R,0.95);
    bwn=XNucleiBinarization(R,bwf,ac);
    
    %% voting based nuclei detection
    d=12;
    Para.rmin=d*0.5;
    Para.rmax=d*1.5;
    Para.bandwidth=d*0.6;
    Para.delta=pi/6;
    Para.Sigma=6;
    [s2]= XVotingClustering(bwn,R,Para);
    
    rs5=round(s2(:,1));cs5=round(s2(:,2));
    %hold on,plot(cs5,rs5,'c+','MarkerSize',13,'LineWidth',2);
    
    
    %% watershed based nuclei segmentation
    fgm4=zeros(size(R));
    ind4=sub2ind(size(R),rs5,cs5);
    fgm4(ind4)=1;
    [bnf,blm] = XWaterShed(bwn,fgm4);
    bnf=bwareaopen(bnf,ac);
    
    %% show the segmentationr results
    B=bwboundaries(bnf);
    hold on,
    for j=1:length(B)
        boundary = B{j};
        boundary(:,1)=smooth(boundary(:,1),3);
        boundary(:,2)=smooth(boundary(:,2),3);
        plot(boundary(:,2), boundary(:,1), 'y-', 'LineWidth', 2);
    end
    
    
end

end
