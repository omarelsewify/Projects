clear, close , clc
V = [0:30];

h = 12.12 - 1.16*V + 11.6*V.^(1/2);

plot(V,h)