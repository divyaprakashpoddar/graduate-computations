clear; clc; close all;

for i = 1:215
    fName = sprintf('U_%d.txt',i*100);
    u = load(fName);
    fName = sprintf('V_%d.txt',i*100);
    v = load(fName);
    
    plot_figures(flipud(u)',flipud(v)','velocity')
    fName = sprintf('newImages/image-%d.png',i);
    saveas(gca,fName)
end

