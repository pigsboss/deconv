function printeps(x,y,Z,x_label,y_label,filename)
close all
figure('Name',filename)
imagesc(x,y,Z)
colorbar('NorthOutside')
axis image
axis xy
xlabel(x_label)
ylabel(y_label)
drawnow
print(filename,'-depsc2')
return