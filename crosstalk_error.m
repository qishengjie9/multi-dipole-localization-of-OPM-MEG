function data_out = crosstalk_error(data_in,sensor_pos,crosstalkError)
n_channels = size(sensor_pos,1);
% add crosstalk error for the MEG data
%crosstalkError =0.05;
if crosstalkError == 0
    data_out = data_in;
else
%% creat the crosstalk matrix
% obtain the distance between the sensors
 distance = zeros(n_channels);
    for i = 1:n_channels
        for j = 1:n_channels
            if(i == j)
                distance(i,j) = 50;
            else
                distance(i,j) = norm(sensor_pos(i,:)-sensor_pos(j,:));
            end        
        end
    end
    distance_min = min(distance,[],'all');
    scale = (distance_min./distance).^3;%set the crosstalk is inversely proportional to the cube of distance
    scale = scale.*crosstalkError;
    scale(logical(eye(size(scale)))) = 1;      % set the diagonal velue of scale to be 1
    tmp = zeros(size(data_in));
% tmp2 = zeros(size(test.data));
% for i = 1:n_channels
%     for j = 1:n_channels
%         for t = 1:size(test.data,2)
%             tmp(i,t) = tmp(i,t)+test.data(j,t)*scale(i,j);
%         end       
%     end 
% end
    for t = 1:size(data_in,2)
        for i = 1:n_channels       
                tmp(i,t) = scale(i,1:n_channels)*data_in(1:n_channels,t);     
        end 
    end
    data_out = tmp;
end

%cross = test.D.data();



