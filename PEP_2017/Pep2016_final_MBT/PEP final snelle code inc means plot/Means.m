h = waitbar(0,'Processing...');
T = 12:2:128;
T(2,:) = T(1,:) - 1;
means = zeros(2,length(T));
waitb = 0;

for P = 1:length(T);    
    waitb = waitb + 1;
    waitbar(waitb/(length(T)),h);
    gridcol = T(1,P);
    gridrow = T(2,P);
    Pinholequality
    means(1,P) = meanplot(1,1);
    means(2,P) = meanplot(1,2);
    means(3,P) = meanplot(1,3);
    means(4,P) = meanplot(1,6);
    means(5,P) = meanplot(1,5);
    means(6,P) = meanplot(1,4);
    disp(means)
    clearvars -except T means waitb h P
end

close(h)

plot(T(1,:),means(:,:),'+');axis([0 140 0 1]);axis square;
title('Mean qualities for different grid sizes')
legend('Mean entire volume','Mean center volume','Median entire volume','Maximum entire volume','Median center volume','Minimum center volume','Location','SouthWest')
xlabel('Number of columns in the grid')
ylabel('Average quality')