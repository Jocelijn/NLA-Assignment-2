nus = [1:10];
G = gallery('grcar',50);

ref_eigenVs = eig(G);


for nu = (1:10)
    hold on;
    for j = (1:10)
        u = rand(50,1);
        v = rand(50,1);
        u = u/norm(u);
        v = v/norm(v);
        SE = G + 10^(-1*nu)*(u*v');
        eigenVs = eig(SE);
        plot(eigenVs,'.k');
    end
    plot(ref_eigenVs,'*r');
    axis([-1,3, -4, 4])
    saveas(gcf,sprintf("EpsilonSpectra/nu_%d_spectra.png",nu));
    clf;
    hold off;
end

    