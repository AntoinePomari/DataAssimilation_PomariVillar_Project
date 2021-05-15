function [mean_error,var_error] = wave1d_accuracy(series_data,observed_data)
error=[];
error=series_data(1:5,:)-observed_data(1:5,2:end);
mean_error=[];
var_error=[];
for i=1:5
    histogram(error,10)
    mean_error(i)=mean(series_data(i,:)-observed_data(i,2:end));
    var_error(i)=var(series_data(i,:)-observed_data(i,2:end));
end
end

