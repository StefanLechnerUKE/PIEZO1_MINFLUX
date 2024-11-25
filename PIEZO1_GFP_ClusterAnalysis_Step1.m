%% Script for the analysis of MINFLUX signals from DNA-Paint labelled GFP-tagged PIEZOs
    % Dependencies: DBSCAN.m; LoadDataFromMFXFile.m; LoadImage2Array.m;
    % SetGFPanalysisParameters.m, PlotRawData.m; ALFAfitSurface
clear all
% close all
%%  %%%%%%%%%%%%%%%%%%%%%%%% – CHOOSE OPTIONS – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    opt.plotRaw             = false; % set to 'true' to plot raw data overview
    opt.PlotSTD             = false;
    opt.FindClusters        = false;
    opt.ShowConfocalImage   = false;    
    aggregateTraces         = true;

ClusterCountNew = 0;
ClusterNameCount = 0;
%% create empty result arrays
    AllCounts = []; AllNumTraces = []; LocPerTraceAll = []; numTracesperHour = [];
StDevHist = [];

%%  %%%%%%%%%%%%%%%%%%%%%%%% – Load data – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      data = readtable('FileList_GFP_OSMO.xlsx');
       % data = readtable('FileList_GFP_CTL_V2.xlsx');

    % for FileIdx = 1:size(data,1);
    for FileIdx = 5 %%  13 is the example in fig 1; idx17 ctl for deep cluster; idx5 for OSMO
        a.myfile = string(data{FileIdx,1});
        a.UpperZ = [data{FileIdx,2}];
        a.LowerZ = [data{FileIdx,3}];
        a.include = [data{FileIdx,4}];
        a.ClusterDistance = 80%0.7*[data{FileIdx,9}];
        a.ClusterMinNumPnts = 100%[data{FileIdx,10}];
        if a.include==0
            continue
        end
        [traces_RAW]=LoadDataFromMFXFile(a.myfile);
    
    %%  %%%%%%%%%%%%%%%%%%%%%%%% – set analysis parameters – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [p]=SetGFPanalysisParameters_V2(); % set analysis paramters

%%  %%%%%%%%%%%%%%%%%%%%%%%% – create empty results matrices – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ClusterNearestNeighborDist = []; NOClusterNearestNeighborDist = []; ALLConfIntensities = [];
    
%%  %%%%%%%%%%%%%%%%%%%%%%%% – plot raw data – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opt.plotRaw
        PlotRawData(traces_RAW,a.myfile);
    end

%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  %%%%%%%%%%%%%%%%%%%%%%%% – filtering – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    traces_FILT = traces_RAW;
    traces_FILT(:,1:3)=1e9*traces_FILT(:,1:3);
    traces_FILT(:,3)=0.7*traces_FILT(:,3);

%% aggregate traces
    if aggregateTraces
    p.loc_per_trace_threshold = 0; traces_Aggregate = []; j=1;
    [w, ~, ~] = unique(traces_FILT(:,4));
    agg = 3;
    for i=1:size(w,1)
    temp = traces_FILT(traces_FILT(:,4)==w(i,1),1:8);
        for k=1:agg:size(temp,1)
            if size(temp,1)<(agg+1)
                continue
            elseif k+(agg-1)>size(temp,1)
                traces_Aggregate(j,1:8) = mean(temp(k:size(temp,1),:),1);
            elseif k+(agg-1)<=size(temp,1)
                traces_Aggregate(j,1:8) = mean(temp(k:k+(agg-1),:));
            end
            j=j+1;
        end
    end
    traces_FILT = traces_Aggregate;
    clear w agg temp;
    end


 % PlotMyData(traces_FILT, Numtraceforeachiteration);
   
 %a.UpperZ = 600e-9;
    %% Z-filtering
        traces_FILT(traces_FILT(:, 3) < 0.7*1e9*a.LowerZ, :)= []; % filter by Z upper limit 160,:) = []; 
        traces_FILT(traces_FILT(:, 3) > 1e9*a.UpperZ, :)= []; % filter by Z lower limit
        % traces_FILT(traces_FILT(:, 3) > 200, :)= []; % filter by Z upper limit 160,:) = []; 
    
    %% XY-filtering
        % % XYcutoffs = [-4975 -2280 11036 8290]; && cutoff für Fig1C
        % XYcutoffs = [-3630 -3380 10640 10400]; %% cutoff for Fig1E
        % XYcutoffs = [-3267 -2929 10581 10251]; %% cutoff für Fig1F
        %  XYcutoffs = [-3080 -2650 8800 8340]; %% cutoff für Fig1F
        % region = 1;
        % traces_FILT(traces_FILT(:, 1) > XYcutoffs(region,2), :)= [];  % optional X-right
        % traces_FILT(traces_FILT(:, 1) < XYcutoffs(region,1), :)= [];  % optional X-left
        % traces_FILT(traces_FILT(:, 2) > XYcutoffs(region,3), :)= [];  % optional Y-upper
        % traces_FILT(traces_FILT(:, 2) < XYcutoffs(region,4), :)= [];  % optional Y-lower
    
    %% efo, cfr filtering
         
        
        traces_FILT = traces_FILT(traces_FILT(:,6) <= p.cfr_threshold, :);

        traces_FILT = traces_FILT(traces_FILT(:,5) <= p.efo_threshold, :); % filter by cfr
        traces_FILT = traces_FILT(traces_FILT(:,5) >= 20000, :);

    %% time filtering 
        traces_FILT = traces_FILT(traces_FILT(:,8) <= p.time_threshold, :);                         % filter by recording time
    
    %% filter by standard deviation
        [uv_tid, ~, id_tid] = unique(traces_FILT(:,4));
        StDevALL=[accumarray(id_tid,traces_FILT(:,1),[],@std) accumarray(id_tid,traces_FILT(:,2),[],@std) accumarray(id_tid,traces_FILT(:,3),[],@std)];
       

        % create arrays for graph coloring purposes only 
        STDforeachiteration = zeros(size(id_tid,1),3);
        for Tidx = 1:size(id_tid,1)
            STDforeachiteration(Tidx,:)=StDevALL(id_tid(Tidx,1),:);
        end
           
        StdMAX = max(STDforeachiteration,[],2); % can be used for loc coloring
        traces_FILT = traces_FILT(StdMAX(:,1)<p.stdev_trace_threshold,:);
        StdMAX = StdMAX(StdMAX(:,1)<p.stdev_trace_threshold,:);
   
        
    %% filter by localisations per trace 
        [uv_tid, ~, id_tid] = unique(traces_FILT(:,4));                                                     % Unique elements and locations in third column  
        n_tid = histcounts(id_tid,size(uv_tid,1)); % How many of each?
        
        % create arrays for graph coloring purposes only 
        LocPerTrace=n_tid';
        Numtraceforeachiteration = id_tid;
        for Tidx = 1:size(id_tid,1)
            Numtraceforeachiteration(Tidx)=LocPerTrace(id_tid(Tidx,1));
        end

        % old filtering code
        traces_FILT = traces_FILT(ismember(traces_FILT(:,4), uv_tid(n_tid > p.loc_per_trace_threshold)),:); % Keep ones with more than threshold.

        % new filtering code
        % traces_filt = traces_filt(Numtraceforeachiteration(:,1)>p.loc_per_trace_threshold,:);
        Numtraceforeachiteration = Numtraceforeachiteration(Numtraceforeachiteration(:,1)>p.loc_per_trace_threshold,:);
        numTracesperHour(FileIdx,1)=size(LocPerTrace,1);
        LocHist = histcounts(LocPerTrace,"BinWidth",1)';
        LocPerTraceAll = cat(1,LocPerTraceAll,LocPerTrace);
        n_tid=n_tid(n_tid(:)>p.loc_per_trace_threshold)';



%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  %%%%% – calculate center of mass for each trace and merge signals from same fluorophore– %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [uv_tid, ~, id_tid] = unique(traces_FILT(:,4)); 
    traces_AVG = [accumarray(id_tid,traces_FILT(:,1),[],@mean) accumarray(id_tid,traces_FILT(:,2),[],@mean) accumarray(id_tid,traces_FILT(:,3),[],@mean)];

    % run DBSCAN to find signals from the same fluorophore – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ClusterIDs = dbscan(traces_AVG(:,1:3),p.MaxDistProtamers,p.MinNumProtamers);
    
    % add column with cluster ID to main coordinates table
    clust_XYZ = [traces_AVG ClusterIDs n_tid];
    
    % split clust_XYZ into 2 submatrices - one with single signals and one with multiple signals per channel 
    singleloc = clust_XYZ(clust_XYZ(:,4) == -1,:);
    multiloc = clust_XYZ(clust_XYZ(:,4) > -1,:);
    
    % calculate averages of clusters (i.e. signals from one channel)
    [uv, ~, id] = unique(multiloc(:,4));
    MyMeanX = [accumarray(id,multiloc(:,1),[],@mean)];
    MyMeanY = [accumarray(id,multiloc(:,2),[],@mean)];
    MyMeanZ = [accumarray(id,multiloc(:,3),[],@mean)];
    ClID = [accumarray(id,multiloc(:,4),[],@mean)];
    numLocs = [accumarray(id,multiloc(:,5),[],@sum)];
    multiloc = [MyMeanX, MyMeanY, MyMeanZ, ClID, numLocs];
    
    % merge submatrices into single matrix
    traces_AVG = cat(1,singleloc(:,1:5),multiloc(:,1:5));


    Xshift = [data{FileIdx,5}];
    Yshift = [data{FileIdx,6}];
    PixelSize = [data{FileIdx,7}];
    myImage = string(data{FileIdx,8});
    myXYZ = 1e-9*traces_FILT;
    [ImageArray, ImageArrayFULL] = LoadImage2Array(myImage,myXYZ,Xshift,Yshift,PixelSize);
 


%%  %%%%%%%%%%%%%%%%%%%%%%%% – identify and highlight clusters - OPTIONAL  – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if opt.FindClusters
            
        %% find and extract clusters from trace RAW data
        ClusterIDsRAW = dbscan(traces_FILT(:,1:3),a.ClusterDistance,a.ClusterMinNumPnts);
         % only consider clusters with more than 12 traces
        for k = 1:max(ClusterIDsRAW)
            
            NumTracesPerCluster = size(unique(traces_FILT(ClusterIDsRAW==k,4)),1);
            if NumTracesPerCluster<10
               ClusterIDsRAW(ClusterIDsRAW(:,1)==k)=-1;
            end
        end

        % create RAW cluster labels
            ClusterIDs_labellist = ClusterIDsRAW(ClusterIDsRAW>-1,:);
            ClusterRAWCoordinates = traces_FILT(ClusterIDsRAW>-1,:);
            ClusterLabelRAWX = [accumarray(ClusterIDs_labellist,ClusterRAWCoordinates(:,1),[],@mean)]+100;
            ClusterLabelRAWY = [accumarray(ClusterIDs_labellist,ClusterRAWCoordinates(:,2),[],@mean)];
            ClusterLabelRAWZ = [accumarray(ClusterIDs_labellist,ClusterRAWCoordinates(:,3),[],@mean)];
            ClusterColorRAWNew = [accumarray(ClusterIDs_labellist,ClusterIDs_labellist(:,1),[],@mean)];

        % write trace RAW loc coordinate data of individual clusters to structure
        for k =1:max(ClusterIDsRAW)
            NumTracesPerCluster = size(unique(traces_FILT(ClusterIDsRAW==k,4)),1);
            if NumTracesPerCluster>=12
                ClusterNameCount=ClusterNameCount+1;
            IndivClusterNames(ClusterNameCount,1) = strcat('FileIdx_',string(FileIdx),'_Clus_',string(k));
       
               IndivClustersRAW.(strcat('FileIdx_',string(FileIdx),'_Clus_',string(k))) = traces_FILT(ClusterIDsRAW==k,:);
            
            %% test code for measuring confocal intensities  
            % ClusterCountNew = ClusterCountNew+1;
            %     tempCluster = traces_FILT(ClusterIDsRAW==k,:);
            %     ClusterDimMinX = 1e-9*min(tempCluster(:,1),[],1);
            %     ClusterDimMaxX = 1e-9*max(tempCluster(:,1),[],1);
            %     ClusterDimMinY = 1e-9*min(tempCluster(:,2),[],1);
            %     ClusterDimMaxY = 1e-9*max(tempCluster(:,2),[],1);
            %     ClusterConfocal = ImageArray(ImageArray(:,1)>ClusterDimMinX-80e-9,:);
            %     ClusterConfocal = ClusterConfocal(ClusterConfocal(:,1)<ClusterDimMaxX+80e-9,:);
            %     ClusterConfocal = ClusterConfocal(ClusterConfocal(:,2)>ClusterDimMinY-80e-9,:);
            %     ClusterConfocal = ClusterConfocal(ClusterConfocal(:,2)<ClusterDimMaxY+80e-9,:);
            %     mean(ClusterConfocal(:,3),1)
            % ClusterConfIntensity(ClusterCountNew,1) = mean(ClusterConfocal(:,3),1);
            % figure;
            % scatter(ClusterConfocal(:,1),ClusterConfocal(:,2),1000,ClusterConfocal(:,3),'square', 'filled')
            % axis equal, axis tight; colormap gray;
            % title(string(k));
            %% end test code
            end
        end
    
        %% find and extract clusters from trace AVG data
        ClusterIDsTraceAVG = dbscan(traces_AVG(:,1:3),p.MaxDistWithinCluster,p.MinNumChannelsPerCluster);
       
        % create distance matrix for nearest neighbor analysis
        DistMatrix=squareform(pdist(traces_AVG(:,1:3)));
        [H, Hidx] = sort(DistMatrix,1);
        NearestNeighbor = H(2,:)';

        % extract results for clusters
        ClusterNumChannels = grouptransform(ClusterIDsTraceAVG,ClusterIDsTraceAVG,@numel);
        ClusterCoordinates = traces_AVG(ClusterIDsTraceAVG>-1,:);
        ClusterNearestNeighborDist = cat(1,ClusterNearestNeighborDist,NearestNeighbor(ClusterIDsTraceAVG>-1,:));
        ClusterIDs_list = ClusterIDsTraceAVG(ClusterIDsTraceAVG>-1,:);
        
        % extract results for signals outside clusters
        NOClusterCoordinates = traces_AVG(ClusterIDsTraceAVG< 0,:);
        NOClusterNearestNeighborDist = cat(1,NOClusterNearestNeighborDist,NearestNeighbor(ClusterIDsTraceAVG< 0,:));
       
        % write trace AVG coordinate data of individual clusters to structure
        for k =1:max(ClusterIDsTraceAVG)
            IndivClustersAVG.(strcat('FileIdx_',string(FileIdx),'_Clus_',string(k))) = traces_AVG(ClusterIDsTraceAVG==k,:);
        end

        % create AVG cluster label
        ClusterLabelX = [accumarray(ClusterIDs_list,ClusterCoordinates(:,1),[],@mean)]+100;
        ClusterLabelY = [accumarray(ClusterIDs_list,ClusterCoordinates(:,2),[],@mean)];
        ClusterLabelZ = [accumarray(ClusterIDs_list,ClusterCoordinates(:,3),[],@mean)];
        ClusterColorNew = [accumarray(ClusterIDs_list,ClusterIDs_list(:,1),[],@mean)];

        figure;
        tiledlayout(1,1);
        % set(gcf,'renderer','Painters');
        % ax = nexttile;
        %     % plot clusters
        %     scatter3(ClusterCoordinates(:,1), ClusterCoordinates(:,2), ClusterCoordinates(:,3),20, ClusterIDs_list(:,1), 'filled'); %'red', 'filled');
        %     colormap(ax,jet);  axis equal; hold on;
        %     % add cluster labels
        %     text(ClusterLabelX(:,1), ClusterLabelY(:,1), ClusterLabelZ(:,1),sprintfc(' %d',ClusterColorNew(:,1)),'FontSize',20)
        %     % plot non-clusters
        %     scatter3(NOClusterCoordinates(:,1), NOClusterCoordinates(:,2), NOClusterCoordinates(:,3),20, [0.6 0.6 0.6], 'filled'); 
       
        ax3 = nexttile;
            scatter3(traces_FILT(:,1), traces_FILT(:,2), traces_FILT(:,3), 10, ClusterIDsRAW(:,1), 'filled');  % color options: ClusterIDs, myXYZ(:,5)
            colormap(ax3,parula); axis equal; title(strcat(a.myfile,' filtered data'));
            text(ClusterLabelRAWX(:,1), ClusterLabelRAWY(:,1), ClusterLabelRAWZ(:,1),sprintfc(' %d',ClusterColorRAWNew(:,1)),'FontSize',20)

    end


%%  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%  %%%%%%%%%%%%%%%%%%%%%%%% – plot processed data – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% %%%%%%%%%%%%%%%%%%%%%%%% - load confocal image and plot overlay with MFX data    

   
    figure;
    tiledlayout(1,2);
    set(gcf,'renderer','Painters');
    ax1 = nexttile;
        scatter(ImageArray(:,1),ImageArray(:,2),300,ImageArray(:,3),'square', 'filled')
        axis equal; axis tight; colormap(ax1,gray); hold on
        scatter(myXYZ(:,1), myXYZ(:,2),10,'red','filled');
    
    ax2 = nexttile;
        % optional code to only visualize the identified clusters
        % traces_filt = traces_filt(ClusterIDsRAW(:,1)>0,:);
        % ClusterIDsRAW = ClusterIDsRAW(ClusterIDsRAW(:,1)>0,:);
        % optional code end
        scatter3(traces_FILT(:,1), traces_FILT(:,2), traces_FILT(:,3),10,traces_FILT(:,3), 'filled');  % color options: ClusterIDs, myXYZ(:,5)
        colormap(ax2,parula); axis equal; title(strcat(a.myfile,' filtered data'));
        % text(ClusterLabelRAWX(:,1), ClusterLabelRAWY(:,1), ClusterLabelRAWZ(:,1),sprintfc(' %d',ClusterColorRAWNew(:,1)),'FontSize',20)
        grid off;
        view(0,90);
if opt.ShowConfocalImage 
    figure;
        scatter(ImageArrayFULL(:,1),ImageArrayFULL(:,2),1000,ImageArrayFULL(:,3),'square', 'filled')
        axis equal
        colormap gray
        hold on
        scatter(myXYZ(:,1), myXYZ(:,2),50,'red','filled');
end

%%  %%%%%%%%%%%%%%%%%%%%%%%% – plot STD dev – %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% StDevHist = []; 
[uv_tid, ~, id_tid] = unique(traces_FILT(:,4));
 StDevALL=[accumarray(id_tid,traces_FILT(:,1),[],@std) accumarray(id_tid,traces_FILT(:,2),[],@std) accumarray(id_tid,traces_FILT(:,3),[],@std)];
       
StDevHist = cat(1, StDevHist, StDevALL(StDevALL(:,1)>0,:));



end % end of main loop

LocHist = histcounts(LocPerTraceAll,"BinWidth",1)';

if opt.PlotSTD
    figure;
    tiledlayout(1,3);
    nexttile
    histogram(StDevHist(:,1),BinWidth=0.5);
    title("STD-X filtered data");

    nexttile
    histogram(StDevHist(:,2),BinWidth=0.5);
    title("STD-Y filtered data");

    nexttile
    histogram(StDevHist(:,3),BinWidth=0.5);
    title("STD-Z filtered data");
end

% ALFAfitSurface(traces_AVG, [],[], 'cubic', 300);

figure;
 % set(gcf,'renderer','Painters');
        scatter3(traces_FILT(:,1), traces_FILT(:,2), traces_FILT(:,3),3,traces_FILT(:,3), 'filled');  % color options: ClusterIDs, myXYZ(:,5)
        colormap(ax2,parula); 
        axis equal; title(strcat(a.myfile,' filtered data'));
        % text(ClusterLabelRAWX(:,1), ClusterLabelRAWY(:,1), ClusterLabelRAWZ(:,1),sprintfc(' %d',ClusterColorRAWNew(:,1)),'FontSize',20)
        grid off;
        % view(0,90);
%         xlim([-3080 -2680])
%         ylim([8340 8740])

% figure;
%  % set(gcf,'renderer','Painters');
%         scatter3(traces_AVG(:,1), traces_AVG(:,2), traces_AVG(:,3),5,traces_AVG(:,6), 'filled');  % color options: ClusterIDs, myXYZ(:,5)
%         colormap(ax2,parula); 
%         axis equal; title(strcat(a.myfile,' filtered data'));
%         % text(ClusterLabelRAWX(:,1), ClusterLabelRAWY(:,1), ClusterLabelRAWZ(:,1),sprintfc(' %d',ClusterColorRAWNew(:,1)),'FontSize',20)
%         grid off;
%         % view(0,90);
% %         xlim([-3080 -2680])
%          zlim([0 1000])

