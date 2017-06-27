function [qualitySas] = Intersections(gridrow,gridcol,voxelCoords,circles,translatedScanning)


% %{
qualitySas = zeros(1,size(voxelCoords,1));
disp(size(voxelCoords,1));
time = zeros(size(voxelCoords,1),1);

for i = 1:size(voxelCoords,1)
    tic
    voxelCoor = voxelCoords(i,1:3);
    singleVoxelPinholes = v2pinhole_scanning_1voxel(translatedScanning, voxelCoor);
    voxelGrid = circle_areaFloSas2(singleVoxelPinholes,gridrow,gridcol);
    %size(voxelGrid)
    %size(circles)

%%Sasha

    voxelIntersect = and(voxelGrid,circles);
    voxelIntersect = any(voxelIntersect,1);
    voxelIntersect = any(voxelIntersect,2);
    qualitySas(i) = sum(voxelIntersect)/size(voxelIntersect,3);
    time(i) = toc;
    ETA_Intersect = (size(voxelCoords,1) - i)*mean(time(1:i))/60;
    if mod(i,10)==0
        disp(ETA_Intersect);
    end
end
% %}
%%Floris

%{
[col,row] = find(voxelGrid);

hit_mattest = zeros(size(row,1),size(circles,3));
[testRow,testCircles] =  meshgrid(row,1:size(circles,3));
[testCol,~] =  meshgrid(col,1:size(circles,3));
testCol = reshape(testCol, [numel(testCol),1]);
testRow =reshape(testRow, [numel(testRow),1]);
testCircles = reshape(testCircles, [numel(testCircles),1]);
index = sub2ind(size(circles),testCol,testRow,testCircles);
hit_mattest = circles(index);

qualtest = sum(hit_mattest(:,:));
qualtest = qualtest > 0;
qualityFlotest = sum(qualtest) / size(circles,3);
%}
%{
%tic
singleVoxelPinholes = v2pinhole_scanning_1voxel(translatedScanning, voxelCoords);
voxelGrid = circle_areaFloSas2(singleVoxelPinholes,gridrow,gridcol);
[row,col] = find(voxelGrid);
hit_mat = zeros(size(col,1),size(circles,3));
for i=1:size(col,1)
    hit_mat(i,:) = circles(row(i),col(i),:);
end
%time-toc

qual = sum(hit_mat(:,:))
qual = qual > 0;
qualityFlo = sum(qual) / size(circles,3);
%}