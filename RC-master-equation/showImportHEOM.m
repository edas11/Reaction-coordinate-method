%clear all;
addpath('RCcore');
addpath('funcScriptData');

densityHEOM = importHEOM('../data/HEOM-dynamics-(IC)11/','(T)300(g-1)100(d)100(l)100(J)100');
t = 1:size(densityHEOM, 1);
dim = size(densityHEOM, 2);

% Figures
for i=1:dim
    for j=1:dim
        subplot(dim,dim,(i-1)*dim+j)
        plot(t, squeeze(real(densityHEOM(:,i,j))), '--b');
        grid on;
        hold on;
        plot(t, squeeze(imag(densityHEOM(:,i,j))), '-r');
    end
end
legend('real','imag');
