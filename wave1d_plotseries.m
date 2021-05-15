function wave1d_plotseries(t,series_data,s,obs_data)
%function wave1d_plotseries(t,series_data,s,obs_data);
%plot timeseries from model and observations
loc_names=s.loc_names;
nseries=length(s.xlocs_waterlevel);
for i=1:nseries
    figure(i);
    grid on
    plot(t,series_data(i,:),'b-');
    hold on
    % plotted in blue our constructed model
    ntimes=min(length(t),size(obs_data,2));
    plot(t(1:ntimes),obs_data(i,2:(ntimes+1)),'r-'); %observations also contain initial time
    title(loc_names{i});
    xlabel('time');
    hold off;
    datetick('x');
    grid on
    legend('Predicted Obs.','Filtered Obs.')
    %print(replace(sprintf('%s.png',loc_names{i}),' ','_'),'-dpng');
end
