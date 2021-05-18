function [mean_error,var_error] = wave1d_accuracy(series_data,observed_data)
error=[];
error=series_data(1:5,:)-observed_data(1:5,2:end);
mean_error=[];
var_error=[];
for ii=1:5
    figure,
    histogram(error,10)
    mean_error(ii)=mean(series_data(ii,:)-observed_data(ii,2:end));
    var_error(ii)=var(series_data(ii,:)-observed_data(ii,2:end));
end
end

