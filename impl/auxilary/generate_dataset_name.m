function [ data_set_name ] = generate_dataset_name( axis_points_name, tensor_size )
% generates a dataset name for grid datasets
    D = numel(tensor_size);
    data_set_name = axis_points_name;
    for d = 1 : D
        data_set_name = [data_set_name sprintf('%d_', tensor_size(d))];
    end
end

