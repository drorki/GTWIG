function  trackersOfCurrentCluster = get_trackers_no(clusterCount, idx, dispFlag)
% Return the names of the trackers that are contained in cluster #clusterCount
uniqueIdx = unique(idx);
clusterName = uniqueIdx(clusterCount);
 trackersOfCurrentCluster = find(idx == clusterName);
 if dispFlag
     disp(['Cluster #' num2str(clusterName) ' contain trackers: ' num2str(trackersOfCurrentCluster(:)') ]);
 end