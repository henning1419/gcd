gx = linspace(0,25,10000)';
gc = zeros(10000,1);

para = [
    1, 0.2, 0.2, 100, 0.5, 0.5;
    2, 0.1, 0.5, 90, 0.5, 0.8;
    4, 0.2, 0.2, 100, 0.5, 0.5;
    4.8, 0.2, 0.2, 100, 0.5, 0.5;
    6, 0.2, 0.2, 100, 0.5, 0.5;
    6.5, 0.2, 0.2, 100, 0.5, 0.5;
    8, 0.2, 0.2, 100, 0.5, 0.5;
    8.3, 0.2, 0.2, 100, 0.5, 0.5;
    10, 0.2, 0.4, 100, 0.5, 0.5;
    11, 0.2, 0.2, 10, 0.5, 0.5;
    14, 0.5, 0.8, 60, 0.5, 0.5;
    15, 0.5, 1.5, 10, 0.5, 0.5;
    ];


for i=1:size(para,1)
    current_para = para(i,:);
    gc = gc + mod_ga(gx,current_para);
    
end
%scale
gc=gc*2;

% noise
abs_noise = 0.1 * randn(size(gc));
gc = gc + abs_noise;
rel_noise = 2 * randn(size(gc)) .* gc./max(max(gc));
gc = gc + rel_noise;


% baseline
base_x = [0 4 10 25];  % [x y] coordinates
base_c = 2*[0 3 2 3];  % [x y] coordinates
baseline = pchip(base_x, base_c, gx);
gc=gc+baseline;


plot(gx,baseline)
hold on
plot(gx,gc)

writematrix([gx gc],[cd filesep 'example_data' filesep 'curve_example.csv'])