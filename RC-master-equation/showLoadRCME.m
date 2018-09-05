% Loads data from .mat file
path = '../data/RCM-fluor-example/';
param = 'RC-(T)300(g-1)100(d)100(l)100(J)100(RCalfa)8(RCdim)8(S)2';
densityLoad = load(strcat( path, param, '.mat' ));
densitySystem = densityLoad.densityReturned;
t = 1:size(densitySystem, 1);
dim = size(densitySystem, 2);

% Figures
for i=1:dim
    for j=1:dim
        subplot(dim,dim,(i-1)*dim+j)
        plot(t, squeeze(real(densitySystem(:,i,j))), '-b');
        grid on;
        hold on;
        plot(t, squeeze(imag(densitySystem(:,i,j))), '-r');
        axis([0 5000 -0.5 0.5]);
    end
end
legend('real','imag');
