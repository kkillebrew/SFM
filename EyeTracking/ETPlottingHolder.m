if options.displayETFigs
    plotNumRows = 5;
    
    % Make x-axis markers for time plots
    % (Used to draw x-axis for each of the 5 plots in
    % the figure)
    clear xAxisMarkers
    maxYVal = 200;
    xAxisMarkers = maxYVal:maxYVal*2:maxYVal*10;
    
    for iI=1:5
        % for iI=1:size(data.et.(options.runType{iB}).etDataOrig,1)
        for iB=1:2   % A/B runs
            
            % Make size values for individual plots
            clear timePlotSize scatterPlotSize
            if iB==1
                timePlotSize(1,:) = [.33 .01 .33 .98];
                scatterPlotSize(1,:) = [0.03 .01 .28 .98];
            elseif iB==2
                for iJ=1:length(data.et.(options.runType{iB}).fixData_diffPosX{iI,iZ})
                    timePlotSize(iJ,:) = [.33 1.01-(.2*iJ) .33 .19];
                    scatterPlotSize(iJ,:) = [.03 1.01-(.2*iJ) .28 .19];
                end
            end
            
            for iZ=1:size(data.et.(options.runType{iB}).etDataOrig,2)   % B/Z days
                % if there is a date for the fix data for this part for this run,
                % then plot the data
                if ~isnan(data.et.(options.runType{iB}).date(iI,iZ))

                    fig1=figure();
                    figSize.baseSize = get(0,'Screensize');   % Base size in pixels
                    figSize.aspectRatio = [10.9849 9.2814];   % Aspect ratio
                    figSize.figSize = [0 0 ...
                        figSize.baseSize(3)...
                        figSize.baseSize(4)];   % Size/postion of fig
                    set(gcf,'color','w')
                    set(gcf,'Position', figSize.figSize)
                    
                    for iJ=1:length(data.et.(options.runType{iB}).fixData_diffPosX{iI,iZ})   % Blocks
                        %% Plot position data in screen coords
                        % subplot(5,3,1+((iJ-1)*3))
                        
                        scatterPlot(iJ) = axes('Parent',fig1);
                        scatter(data.et.(options.runType{iB}).fixData_avePosX{iI,iZ}{iJ},...
                            data.et.(options.runType{iB}).fixData_avePosY{iI,iZ}{iJ});
                        hold on
                        plot([options.screenSize(3)/2 options.screenSize(3)/2],...
                            [options.screenSize(4)/2 options.screenSize(4)/2],'.r');
                        text(options.screenSize(3)*.05,options.screenSize(4)*.4,...
                            sprintf('%s\n%.2f','X Correction Used:',...
                            data.et.(options.runType{iB}).aveDistFromFix_PosX{iI,iZ}(iJ)));
                        text(options.screenSize(3)*.05,options.screenSize(4)*.15,...
                            sprintf('%s\n%.2f','Y Correction Used:',...
                            data.et.(options.runType{iB}).aveDistFromFix_PosY{iI,iZ}(iJ)));
                        xlim([options.screenSize(1) options.screenSize(3)])
                        ylim([options.screenSize(2) options.screenSize(4)])
                        hold off
                        
                        box off
%                         title(['Fixation points (block: ', ...
%                             num2str(iJ),')'],...
%                             'fontsize',10)
                        
                        set(scatterPlot(iJ),'xtick',[],...
                            'ytick',[],...
                            'Position',scatterPlotSize(iJ,:))
                        ylabel(['Block: ', num2str(iJ)])
                        xlabel([])
                    end
                        
                    %% Plot the heat map
                    % subplot(5,3,[3,6])
                    heatPlot(1) = axes('Parent',fig1);
                    % Average data across blocks
                    imagesc(nanmean(cat(3,data.et.(options.runType{iB}).fixData_heatMap{iI,iZ}{:}),3)');
                    hold on
                    plot([options.screenSize(3)/2 options.screenSize(3)/2],...
                        [options.screenSize(4)/2 options.screenSize(4)/2],'.r');
                    hold off
                    set(gca,'YDir','normal')
                    text(options.screenSize(3)*.1,options.screenSize(4)*.25,...
                            sprintf('%s\n%.2f','Ave Euc. Distance From Fix',...
                            nanmean(data.et.B.stats.aveEucDisFromFix{iI,iZ})));
                    text(options.screenSize(3)*.1,options.screenSize(4)*.1,...
                        sprintf('%s\n%d','Num Valid Fixations: ',...
                        data.et.(options.runType{iB}).nonNanCounter{iI,iZ}(iJ)))
                    title(['Proportion of eye gaze position' ...
                        sprintf('\n%s%d%s%d%s',' (P',data.et.(options.runType{iB}).subjID(iI),', ',iZ,' run)')],...
                        'fontsize',10)
                    set(heatPlot(1),'xtick',[],...
                        'ytick',[],...
                        'OuterPosition',[.66 .5 .33 .5])
                    
                    % subplot(5,3,[12,15])
                    heatPlot(2) = axes('Parent',fig1);
                    % Average data across blocks
                    imagesc(nanmean(cat(3,data.et.(options.runType{iB}).corFixData_heatMap{iI,iZ}{:}),3)');
                    hold on
                    plot([options.screenSize(3)/2 options.screenSize(3)/2],...
                        [options.screenSize(4)/2 options.screenSize(4)/2],'.r');
                    hold off
                    set(gca,'YDir','normal')
                    text(options.screenSize(3)*.1,options.screenSize(4)*.25,...
                            sprintf('%s\n%.2f','Ave Euc. Distance From Fix (cor.)',...
                            nanmean(data.et.B.stats.corrAveEucDisFromFix{iI,iZ})));
                    text(options.screenSize(3)*.1,options.screenSize(4)*.1,...
                        sprintf('%s\n%d','Num Valid Fixations: ',...
                        data.et.(options.runType{iB}).nonNanCounterCorr{iI,iZ}(iJ)))
                    set(gca,'XColor','k','YColor','k')
                    title(['Proportion of eye gaze position (corr.)' ...
                        sprintf('\n%s%d%s%d%s',' (P',data.et.(options.runType{iB}).subjID(iI),', ',iZ,' run)')],...
                        'fontsize',10)
                    set(heatPlot(2),'xtick',[],...
                        'ytick',[],...
                        'OuterPosition',[.66 0 .33 .5])
                        
                        
                        
                    %% Plot relative eye position in lx, rx, ly, ry
                    % Find length and slpit evenly into 5 rows to plot
                    for iJ=1:length(data.et.(options.runType{iB}).fixData_diffPosX{iI,iZ})   % Blocks
                        
                        % Make x-axis
                        % (Will be same x vals for all plots)
                        xAxisLims = linspace(1,length(data.et.(options.runType{iB}).fixData_diffPosX{iI,iZ}{iJ}),plotNumRows+1);
                        
                        % Find x-axis labels based on time
                        xAxisTime = linspace(1,data.et.(options.runType{iB}).fixDataTotalTime{iI,iZ}(iJ),plotNumRows+1);
                        
                        % Plot 5 sets of lines w/ full eye position for
                        % each block
                        timePlot(iJ) = axes('Parent',fig1);
                        for iK=1:plotNumRows
                            
                            % Center
                            h1=plot([xAxisLims(1) xAxisLims(2)],[xAxisMarkers(iK) xAxisMarkers(iK)],'k');
                            h1.Annotation.LegendInformation.IconDisplayStyle = 'off';
                            hold on
                            % X Right eye
                            plot(xAxisLims(1):xAxisLims(2),...
                                data.et.(options.runType{iB}).fixData_diffPosX{iI,iZ}{iJ}(...
                                xAxisLims(iK):xAxisLims(iK+1),1)+xAxisMarkers(iK),...
                                'Color',[0 .9 0]);
                            % Y Right eye
                            h2 = plot(xAxisLims(1):xAxisLims(2),...
                                data.et.(options.runType{iB}).fixData_diffPosY{iI,iZ}{iJ}(...
                                xAxisLims(iK):xAxisLims(iK+1),1)+xAxisMarkers(iK),...
                                'Color',[0 .6 0]);
                            h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
                            % X Left eye
                            plot(xAxisLims(1):xAxisLims(2),...
                                data.et.(options.runType{iB}).fixData_diffPosX{iI,iZ}{iJ}(...
                                xAxisLims(iK):xAxisLims(iK+1),2)+xAxisMarkers(iK),...
                                'Color',[0 0 .9])
                            % Y Left eye
                            plot(xAxisLims(1):xAxisLims(2),...
                                data.et.(options.runType{iB}).fixData_diffPosY{iI,iZ}{iJ}(...
                                xAxisLims(iK):xAxisLims(iK+1),2)+xAxisMarkers(iK),...
                                'Color',[0 0 .6])
                        
                            
                            %% Plot blink marker
                            % Left eye blink
                            % If there are left blinks for this trial
                            if isfield(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ},'left')
                                % If there are blinks for this trial
                                if ~isempty(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left)
                                    for iL = 1:length(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left)
                                        % If the blink time falls w/in the plot time limit plot it
                                        if data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left(iL) > xAxisTime(iK) &&...
                                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left(iL) < xAxisTime(iK+1)
                                            
                                            % Scale x values (eyeblink times) based on xAxisTime(iK) to ensure each
                                            % point is plotted relative to x axis '0'. Scale y values relative to 
                                            % individual x axis markers.
                                            patch([data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.left(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.left(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.left(iL)/2]-...
                                                xAxisTime(iK)/2,...
                                                [maxYVal maxYVal -maxYVal -maxYVal]+xAxisMarkers(iK),...
                                                [0 1 0],...
                                                'FaceAlpha',.5,...
                                                'EdgeColor','none')
                                        end
                                    end
                                end
                            end
                            % Right eye blink
                            % If there are rightblinks for this trial
                            if isfield(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ},'right')
                                % If there are blinks for this trial
                                if ~isempty(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right)
                                    for iL = 1:length(data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right)
                                        % If the blink time falls w/in the plot time limit plot it
                                        if data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right(iL) > xAxisTime(iK) &&...
                                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right(iL) < xAxisTime(iK+1)
                                            patch([data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.right(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkEndTime{iI,iZ}{iJ}.right(iL)/2 ...
                                                data.et.(options.runType{iB}).blinkStartTime{iI,iZ}{iJ}.right(iL)/2]-...
                                                xAxisTime(iK)/2,...
                                                [maxYVal maxYVal -maxYVal -maxYVal]+xAxisMarkers(iK),...
                                                [0 0 1],...
                                                'FaceAlpha',.5,...
                                                'EdgeColor','none')
                                        end
                                    end
                                end
                            end
                            
%                             %% Plot sacc marker
%                             % Left eye sacc
%                             % If there are left saccs for this trial
%                             if isfield(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ},'left')
%                                 % If there are blinks for this trial
%                                 if ~isempty(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left)
%                                     for iL = 1:length(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left)
%                                         % If the sacc time falls w/in the plot time limit plot it
%                                         if data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left(iL) > xAxisTime(iK) &&...
%                                                 data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left(iL) < xAxisTime(iK+1)
%                                             patch([data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left(iL)/2 ...
%                                                 data.et.(options.runType{iB}).saccEndTime{iI,iZ}{iJ}.left(iL)/2 ...
%                                                 data.et.(options.runType{iB}).saccEndTime{iI,iZ}{iJ}.left(iL)/2 ...
%                                                 data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.left(iL)/2]-...
%                                                 xAxisTime(iK)/2,...
%                                                 [maxYVal maxYVal -maxYVal -maxYVal]+xAxisMarkers(iK),...
%                                                 [0 1 0],...
%                                                 'FaceAlpha',.1,...
%                                                 'EdgeColor','none')
%                                         end
%                                     end
%                                 end
%                             end
%                             % Right sacc
%                             % If there are right saccs for this trial
%                             if isfield(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ},'right')
%                                 % If there are saccs for this trial
%                                 if ~isempty(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right)
%                                     for iL = 1:length(data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right)
%                                         % If the sacc time falls w/in the plot time limit plot it
%                                         if data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right(iL) > xAxisTime(iK) &&...
%                                                 data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right(iL) < xAxisTime(iK+1)
%                                             patch([data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right(iL)/2 ...
%                                                 data.et.(options.runType{iB}).saccEndTime{iI,iZ}{iJ}.right(iL)/2 ...
%                                                 data.et.(options.runType{iB}).saccEndTime{iI,iZ}{iJ}.right(iL)/2 ...
%                                                 data.et.(options.runType{iB}).saccStartTime{iI,iZ}{iJ}.right(iL)/2]-...
%                                                 xAxisTime(iK)/2,...
%                                                 [maxYVal maxYVal -maxYVal -maxYVal]+xAxisMarkers(iK),...
%                                                 [0 0 1],...
%                                                 'FaceAlpha',.1,...
%                                                 'EdgeColor','none')
%                                         end
%                                     end
%                                 end
%                             end
%                             
                            %% Plot behavioral events
                            % Actual switch times
                            for iL = 1:length(data.behav.controlSwitchTimes)
                                % If the event time falls w/in the plot time limit plot it
                                if data.behav.controlSwitchTimes(iL) > xAxisTime(iK) &&...
                                        data.behav.controlSwitchTimes(iL) < xAxisTime(iK+1)
                                    plot([data.behav.controlSwitchInd(iL)...
                                        data.behav.controlSwitchInd(iL)]-...
                                                xAxisTime(iK)/2,...
                                        [maxYVal -maxYVal]+xAxisMarkers(iK),...
                                        'Color',[0 0 0])
                                end
                            end
                            
                            % Responses
                            % If there are events for this trial
                            if ~isempty(data.behav.(options.runType{iB}).flipTimes{iI,iZ}{iJ})
                                for iL = 1:length(data.behav.(options.runType{iB}).flipTimes{iI,iZ}{iJ})
                                    % If the event time falls w/in the plot time limit plot it
                                    if data.behav.(options.runType{iB}).flipTimes{iI,iZ}{iJ}(iL) > xAxisTime(iK) &&...
                                            data.behav.(options.runType{iB}).flipTimes{iI,iZ}{iJ}(iL) < xAxisTime(iK+1)
                                        plot([data.behav.(options.runType{iB}).flipTimesIdx{iI,iZ}{iJ}(iL)...
                                            data.behav.(options.runType{iB}).flipTimesIdx{iI,iZ}{iJ}(iL)]-...
                                                xAxisTime(iK)/2,...
                                            [maxYVal -maxYVal]+xAxisMarkers(iK),...
                                            'Color',[1 0 0])
                                    end
                                end
                            end
                            
                            xlim([xAxisLims(1) xAxisLims(2)])
                            ylim([0 xAxisMarkers(end)+maxYVal]);
                            set(timePlot(iJ),'XColor','k','YColor','k')
                        end
                        set(timePlot(iJ),'xtick',[],...
                            'ytick',[],...
                            'Position',timePlotSize(iJ,:))
                        if iJ==1
                            legend({'Right','Left'},'Location','NorthWest')
                        end
                    end
                end
            end
        end
    end
end