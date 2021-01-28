function wykres(X,Y,a)
% X - wylosowane punkty (w 3D)
% Y - etykiety punktów
% a - wspó³czynniki liniowego separatora binarnego (4 wsp dla p³aszcyzny rozcinaj¹cej)
%     a1*x1+a2*x2+a3*x3+a4=0

n=size(X,1);    % n=3 punkty w R3
m=size(X,2);    % m=100 punktów treningowych

if n==3 && a(3)~=0
    grid on
    hold on

    for j=1:m
        if Y(j)>0
            plot3(X(1,j),X(2,j),X(3,j),'r.');
        else
            plot3(X(1,j),X(2,j),X(3,j),'b.');
        end
    end


    x=-10:0.1:10;
    [X1,X2]=meshgrid(x);

    X3=-(a(1).*X1+a(2).*X2+a(4))/a(3);
    surf(X1,X2,X3);

end

end