function printeps(x,y,Z,x_label,y_label,filename)
figure('Name',filename)
imagesc(x,y,Z)
truesize
axis equal
axis xy
xlabel(x_label)
ylabel(y_label)
drawnow
print(filename,'-depsc2')
return