function plotTankHeights(T, H, title_)

fig = figure;

subplot(2,2,3);
plot(T, H(:,1));
title('Tank 1')

subplot(2,2,4);
plot(T, H(:,2));
title('Tank 2')

subplot(2,2,1);
plot(T, H(:,3));
title('Tank 3')

subplot(2,2,2);
plot(T, H(:,4));
title('Tank 4')

sgtitle(fig, title_)

end