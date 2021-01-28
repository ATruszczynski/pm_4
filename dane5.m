function [X,Y] = dane5

% Losowe dane z kostki 3d ze współrzędnymi w zakresie (-10, 10)

rng(1111)

a = 10;
n = 1000;

X = [];
Y = [];

for i = 1:n
    point = [];
    for j = 1:3
        point = [point, (rand()-0.5)*2*a];
    end
    point = point.';
    X = [X, point];
    
    col = randi(2);
    if col == 2
        col = -1;
    end
    
    Y = [Y, col];
end

end