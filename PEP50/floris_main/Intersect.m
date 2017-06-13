function I = Intersect(sin_grid, pin_grid)
    %the purpose of the function is to check if the two submitted matrices
    %(functions on the grids), intersect or not.
    [col,row] = find(pin_grid);
    hit_mat = zeros(size(row,1),size(sin_grid,3));
    for i=1:size(row,1)
        hit_mat(i,:) = sin_grid(col(i),row(i),:);
    end
    qual = sum(hit_mat(:,:));
    qual = qual > 0;
    quality = sum(qual) / size(sin_grid,3)
end