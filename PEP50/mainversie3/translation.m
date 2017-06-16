function pinholes = translation(scanning_input, positions, gridrow, gridcol,continu)

stops = size(positions,1);
if continu == true
    numberSteps = max(gridrow,gridcol);
else 
    numberSteps = 1;
end
steps = linspace(0,1,numberSteps);
steps = repmat(steps',[1,3]);
steps = repmat(steps,[stops,1]);
direction = diff(positions,1);
direction = repmat(direction,[numberSteps*stops,1]);

steppedTranslation = direction.*steps;
steppedPositions = steppedTranslation+repmat(positions,[numberSteps, 1]);
steppedPositions = repmat(steppedPositions,[size(scanning_input,1),1]);
pinholes = repmat(scanning_input,[numberSteps*stops,1]);
pinholes(:,1:3) = pinholes(:,1:3)+steppedPositions;

