
% These are the UY equiv WL data for fixed mosaics (5.4 center).
% Entered from the .txt files by hand by DHB, 9/30/24 and 10/1/24

LMRatio(1) = 0.67;
UYDataVSize{1} = [ ...
    1.5	616	NaN
    2.5	583	597
    3.4	570	585
    4.4	566	580
    5.4	569	581
    7.4	569	582
    10.3 570 582 ...
    ];

figure; clf; hold on;
plot(UYDataVSize(:,1),UYDataVSize(:,2),'ro','MarkerFaceColor','r','MarkerSize',12);
plot(UYDataVSize(:,1),UYDataVSize(:,3),'bo','MarkerFaceColor','b','MarkerSize',12);
plot(UYDataVSize(:,1),UYDataVSize(:,2),'r','LineWidth',2);
plot(UYDataVSize(:,1),UYDataVSize(:,3),'b','LineWidth',2);
xlabel('Square size (min)','FontSize',18);
ylabel('Unique Yellow Equiv Wl (nm)','FontSize',18);
title('Unique yellow for 2:1 L:M mosaic','FontSize',20);
legend({'570 nm reference yellow','580 nm reference yellow'},'FontSize',14);