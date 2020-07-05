function [filtered_params] = fun_plot_fit_data(final_params, plotAllEvolutionMat, bigChallengeMat, bigReChallengeMat, plotLowAllEvolutionMat, bigLowChallengeMat, bigReLowChallengeMat,...
                                            primary_community, secondary_community, cytokine_formatted, total_time_in_hours, scatterMat, meta_data_for_file_name, search, filterValue)

    font = 'Times New Roman';                                    
    if (search==1)
        csvwrite(strcat('final_params_',meta_data_for_file_name,'.csv'),final_params)  
    end   
    
    colorMe='c';
    column_names={'alpha', 'beta', 'gamma', 'gamma_2', 'beta_2', 'rmsd'};
    %Filter results based on deviation score
    ind1 = final_params(:,6)<filterValue; %grab indices along third column where condition
    %keep rows with least variance and send to plotting function
    bestFitParams = final_params(ind1,:);
    %1000/1000
    bestFitplotAllEvolutionMat = plotAllEvolutionMat(ind1,:);
    bestFitbigChallengeMat = bigChallengeMat(ind1,:);
    bestFitbigReChallengeMat = bigReChallengeMat(ind1,:);
    %1000/1000 pies
    bestFitprimary_community = primary_community(ind1,:);
    bestFitsecondary_community = secondary_community(ind1,:);
    %10/1000
    bestFitplotLowAllEvolutionMat = plotLowAllEvolutionMat(ind1,:);
    bestFitbigLowChallengeMat = bigLowChallengeMat(ind1,:);
    bestFitbigReLowChallengeMat = bigReLowChallengeMat(ind1,:);

     [result_rows, columns]=size(bestFitParams);

     if result_rows > 0

            plotAll=figure('DefaultTextFontName', font, 'DefaultAxesFontName', font);
            set(gca,'fontname','Times New Roman')
            subplot(2,1,1)
            plot(total_time_in_hours, bestFitplotAllEvolutionMat, colorMe)
            xlabel('Time (in hr)');
            %ylim([0 1]);
            hold on

            %Descriptive stats on all matrices
            bigChallengeMatAvg = mean(bestFitbigChallengeMat);
            bigChallengeMatSD = 1.96*std(bestFitbigChallengeMat);
            bigReChallengeMatAvg = mean(bestFitbigReChallengeMat);
            bigReChallengeMatSD = 1.96*std(bestFitbigReChallengeMat);

            bigMatAvg=[bigChallengeMatAvg bigReChallengeMatAvg];
            bigMatSD=[bigChallengeMatSD bigReChallengeMatSD];
            e=errorbar(0:1:47,bigMatAvg,bigMatSD);
            e.Color='k';

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%Mark experimental data%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %     e2loTol8, e2loTol12, e2loTol16
            %     e2Tol8, e2Tol12, e2Tol16
            %     e2chall8, e2chall12, e2chall16

            scatter(4,scatterMat(1,1),60,'+','MarkerEdgeColor','r',...
                          'MarkerFaceColor','r',...
                          'LineWidth',1.5)
            scatter(8,scatterMat(1,2),60,'+','MarkerEdgeColor','r',...
                          'MarkerFaceColor','r',...
                          'LineWidth',1.5)
            scatter(12,scatterMat(1,3),60,'+','MarkerEdgeColor','r',...
                          'MarkerFaceColor','r',...
                          'LineWidth',1.5)
            scatter(28,scatterMat(3,1),60,'*','MarkerEdgeColor','r',...
                  'MarkerFaceColor','r',...
                  'LineWidth',1.5)
            scatter(32,scatterMat(3,2),60,'*','MarkerEdgeColor','r',...
                  'MarkerFaceColor','r',...
                  'LineWidth',1.5)
            scatter(36,scatterMat(3,3),60,'*','MarkerEdgeColor','r',...
                  'MarkerFaceColor','r',...
                  'LineWidth',1.5)

             hold off

            %figure(2);
            subplot(2,1,2)
            plot(total_time_in_hours, bestFitplotLowAllEvolutionMat, colorMe)
            xlabel('Time (in hr)');
            %ylim([0 1]);
            hold on

            %Descriptive stats on all matrices
            bigLowChallengeMatAvg = mean(bestFitbigLowChallengeMat);
            bigLowChallengeMatSD = 1.96*std(bestFitbigLowChallengeMat);
            bigReLowChallengeMatAvg = mean(bestFitbigReLowChallengeMat);
            bigReLowChallengeMatSD = 1.96*std(bestFitbigReLowChallengeMat);

            bigLowMatAvg=[bigLowChallengeMatAvg bigReLowChallengeMatAvg];
            bigLowMatSD=[bigLowChallengeMatSD bigReLowChallengeMatSD];
            e=errorbar(0:1:47,bigLowMatAvg,bigLowMatSD);
            e.Color='k';

            scatter(28,scatterMat(2,1),60,'*','MarkerEdgeColor','r',...
                  'MarkerFaceColor','r',...
                  'LineWidth',1.5)
            scatter(32,scatterMat(2,2),60,'*','MarkerEdgeColor','r',...
                  'MarkerFaceColor','r',...
                  'LineWidth',1.5)
            scatter(36,scatterMat(2,3),60,'*','MarkerEdgeColor','r',...
                  'MarkerFaceColor','r',...
                  'LineWidth',1.5)

             hold off
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%Don't buy, make pies%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            communityPie=figure('DefaultTextFontName', font, 'DefaultAxesFontName', font)
            set(gca,'fontname','Times New Roman')
            %subplot(2,2,3)
            %RGB colour definitions. 
            posCol = [171 16 16];
            nrsCol = [41 114 70];
            nrpsCol = [229 186 69];
            %Divide by 255 for matlab colour codes
            colormap([posCol/255;      %//maroon
                     0 0 0;      %// black 
                      nrsCol/255;      %// green/nrs
                      nrpsCol/255]);
            title1Text=strcat('Post primary 16hr', cytokine_formatted);   
            title2Text=strcat('Post secondary', cytokine_formatted);
            labels={'Positive: '; 'Negative: '; 'NRS: '; 'NRPS: '};   

            subplot(2,1,1)
            explode = [0 1 0 1 ];
            pieData=mean(bestFitprimary_community);
            pieData(abs(pieData)==0)=1e-3;
            p = pie(pieData, explode);
            pText = findobj(p,'Type','text');
            percentValues = get(pText,'String');
            combinedtxt = strcat(labels,percentValues);
            pText(1).String = combinedtxt(1);
            pText(2).String = combinedtxt(2);
            pText(3).String = combinedtxt(3);
            pText(4).String = combinedtxt(4);
            % ax = gca();
            % set(ax,'fontsize', 24);
            %        pos = get(gca, 'Position');
            %     pos(1) = 0.055;
            %     pos(3) = 0.9;
            %     set(gca, 'Position', pos)
            title(title1Text)

            subplot(2,1,2)
            explode = [0 1 0 1 ];
            %take first 24 items
            pieData=mean(bestFitsecondary_community);
            pieData(abs(pieData)==0)=1e-3;
            %pieData=prot_re_lps_timeCommunity(:,end);
            p = pie(pieData, explode);
            pText = findobj(p,'Type','text');
            percentValues = get(pText,'String');
            combinedtxt = strcat(labels,percentValues);
            pText(1).String = combinedtxt(1);
            pText(2).String = combinedtxt(2);
            pText(3).String = combinedtxt(3);
            pText(4).String = combinedtxt(4);
            % ax = gca();
            % set(ax,'fontsize', 24);
            %  pos = get(gca, 'Position');
            %     pos(1) = 0.055;
            %     pos(3) = 0.9;
            %     set(gca, 'Position', pos)
            title(title2Text)
            %legend(labels)
            filtered_par=figure('DefaultTextFontName', font, 'DefaultAxesFontName', font)
            boxplot(bestFitParams, column_names)
            %set(0,'DefaultFigureVisible','off');
            %title(strcat('fiiltered_params_',meta_data_for_file_name))
             % print all the jazz
            path2save='D:\GoogleDrive\silence\PhD\MATLAB\8_state_models\May_2019\Individual_model_2_non_responsive_states\results\';
            set(filtered_par,'PaperPositionMode','auto'); 
            set(filtered_par,'PaperSize',[25 12]); %set the paper size to what you want  
            print(filtered_par, strcat(path2save,'fiiltered_params_',meta_data_for_file_name),'-dpdf') % then print it

            set(communityPie,'PaperPositionMode','auto'); 
            set(communityPie,'PaperSize',[25 12]); %set the paper size to what you want  
            print(communityPie,strcat(path2save,'communityPie_',meta_data_for_file_name),'-dpdf') % then print it

            set(plotAll,'PaperPositionMode','auto'); 
            set(plotAll,'PaperSize',[25 12]); %set the paper size to what you want  
            print(plotAll,strcat(path2save,'plotAll_',meta_data_for_file_name),'-dpdf') % then print it
     end     
     filtered_params = bestFitParams;
end
 