function [mean_error,var_error,cov_error,quant25,median,quant75] = wave1d_accuracy(series_data,observed_data)
error=[];
error=series_data(1:5,:)-observed_data(1:5,2:end);
mean_error=[];
var_error=[];
cov_error=[];
names={'Cadzand','Vlissingen','Terneuzen','Hansweert','Bath'};
cov_error=cov(error');
for ii=1:5
    figure
    plot(error(ii,:))
    title(names{ii})
    xlabel('time');
    ylabel('Error');
    grid on
    mean_error(ii)=mean(series_data(ii,:)-observed_data(ii,2:end));
    var_error(ii)=var(series_data(ii,:)-observed_data(ii,2:end));
    a=sort(error(ii,:));
    quant25(ii)=a(floor(length(a)/4));
    median(ii)=a(floor(length(a)/2));
    quant75(ii)=a(floor(length(a)*3/4));
end
end

