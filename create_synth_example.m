gx = linspace(0,40,20000)';
gc = zeros(20000,1);

para = [
    16, 0.1, 0.1, 1000, 0.5, 0.5;
    16.6, 0.05, 0.2, 900, 0.5, 0.8;
    19, 0.1, 0.1, 1000, 0.5, 0.5;
    19.3, 0.1, 0.1, 1000, 0.5, 0.5;
    21, 0.1, 0.1, 1000, 0.5, 0.5;
    21.2, 0.1, 0.1, 1000, 0.5, 0.5;
    23, 0.1, 0.1, 1000, 0.5, 0.5;
    23.15, 0.1, 0.1, 1000, 0.5, 0.5;
    25, 0.1, 0.15, 1000, 0.5, 0.5;
    25.5, 0.1, 0.1, 200, 0.5, 0.5;
    29, 0.2, 0.4, 600, 0.5, 0.5;
    29.6, 0.2, 1, 100, 0.5, 0.5;
    35, 0.1, 0.1, 1000, 0.5, 0.5;
    35.3, 0.1, 0.1, 1000, 0.5, 0.5;
    35.7, 0.12, 0.12, 1200, 0.5, 0.5;
    36, 0.1, 0.1, 800, 0.5, 0.5;
    36.4, 0.15, 0.15, 1000, 0.5, 0.5;
    ];


for i=1:size(para,1)
    current_para = para(i,:);
    gc = gc + mod_ga(gx,current_para);
    
end

% noise
abs_noise = 2 * randn(size(gc));
gc = gc + abs_noise;
rel_noise = 2 * randn(size(gc)) .* gc./max(max(gc));
gc = gc + rel_noise;


% baseline
base_x = [0 10 25 30 40];  % [x y] coordinates
base_c = 30*[2 5 0.2 0.5 3];  % [x y] coordinates
baseline = pchip(base_x, base_c, gx);
gc=gc+baseline;


plot(gx,baseline)
hold on
plot(gx,gc)

writematrix([gx gc],[cd filesep 'example_data' filesep 'curve_example.csv'])