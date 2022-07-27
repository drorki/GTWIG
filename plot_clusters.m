function plot_clusters(GT, idx, figName, markRectangleFlag)
% Input:
% 1) GT - array of whistle traces.
% 2) idx - cluster index for each whistle trace
% 3) figName - name for figure. You can also use a valid figure handle instead to plot on an existing figure
% 4) markRectangleFlag - 1/0 to mark/don't mark clusters by a bounding rectangle
%

markTrackerNo = 1; % 1/0 to mark/don't mark traces by theis number
useDistinctLineStyle = 1; % 1/0 to use/don't use multiple line styles

%% --- Plot clusters ----
colorMat = [0  0.4470  0.7410; ...
    0.8500  0.3250  0.0980; ...
    0.9290  0.6940  0.1250; ...
    0.4940  0.1840  0.5560; ...
    0.4660  0.6740  0.1880; ...
    0.3010  0.7450  0.9330; ...
    0.6350  0.0780  0.1840];

lineStyleCell = {'-','--','-.',':'};

if ishandle(figName) 
    figure(figName);   hold on % Use existing figure
else
    figure('Name', figName);  % Create a new figure
end

nClusters = numel(unique(idx));
colorInd = 0;
lineStyleInd = 0;

for clusterCount = 1:nClusters
    trackersOfCurrentCluster = get_trackers_no(clusterCount, idx, 1);
    if ~isempty(trackersOfCurrentCluster)

        colorInd = colorInd + 1;
        if colorInd > size(colorMat, 1)
            colorInd = 1;
        end
        lineStyleInd =lineStyleInd  +1;
        if lineStyleInd > length(lineStyleCell)
            lineStyleInd = 1;
        end
        tMin = 1e12; tMax = 0;
        fMin = 1e12; fMax = 0;
        
        
        for kk =1:length(trackersOfCurrentCluster)
            
            if useDistinctLineStyle
                plot(GT(trackersOfCurrentCluster(kk)).time, GT(trackersOfCurrentCluster(kk)).freq/1000, ...
                    'LineStyle', lineStyleCell{lineStyleInd}, 'Color', colorMat( colorInd, :) ,'LineWidth',2); hold on
            else
                plot(GT(trackersOfCurrentCluster(kk)).time, GT(trackersOfCurrentCluster(kk)).freq/1000, '.-', ...
                    'Color', colorMat( colorInd, :) ); hold on
            end
            if markTrackerNo
                text(mean( GT(trackersOfCurrentCluster(kk)).time ), mean(GT(trackersOfCurrentCluster(kk)).freq/1000), ...
                    num2str(trackersOfCurrentCluster(kk)) )
            end
            % Parameter for bounding box
            tMin = min(tMin, GT(trackersOfCurrentCluster(kk)).time(1));
            tMax = max(tMax, GT(trackersOfCurrentCluster(kk)).time(end));
            fMin =  min(fMin, min(GT(trackersOfCurrentCluster(kk)).freq/1000));
            fMax = max(fMax, max(GT(trackersOfCurrentCluster(kk)).freq/1000));
        end
        if numel(trackersOfCurrentCluster)>1 && markRectangleFlag
            rectangle('Position',[tMin   fMin   tMax-tMin  fMax-fMin], 'LineStyle', '-', 'EdgeColor', 'r')
        end
    end
end
xlabel('time [s]'); ylabel('frequency [kHz]');
nTrackers = length(GT);
title([num2str(nTrackers) ' Original Traces \rightarrow ' num2str(nClusters) ' Segmented Traces'])



