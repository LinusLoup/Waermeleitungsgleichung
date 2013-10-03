function waermeleitungsgleichung

%% Konfiguration
%border = 4;
x_min = -6;
x_max = 6;
y_min = -4;
y_max = 4;
h = 0.1;
dt = 0.001;
%c = 1.25;

[x,y] = meshgrid(x_min:h:x_max,y_min:h:y_max);


%% Startwerte

[mx,nx] = size(x);
[my,ny] = size(y);

inner_x = x(2:mx-1,2:nx-1);
inner_y = y(2:my-1,2:ny-1);

u = zeros(size(x));
U = u;                          % Zwischenlösung

index_x = 1:nx;
index_y = 1:mx;

center_x = length(index_x)/2;
center_y = length(index_y)/2;

    % Temperaturleitungskoeffizient-Matrix für einen kompletten Indexvektor
    function z = temp(ind_y,ind_x)
        
        z = zeros(length(ind_y),length(ind_x));
        for n = 1:length(ind_x)
            for m = 1:length(ind_y)
                if (abs(center_y - m) < 5) && (abs(center_x - n) < 5)
                    z(m,n) = 98.8e-2;
                else
                    z(m,n) = 20e-2;
                end
            end
        end
    end

% Temperaturleitungskoeffizient für das Gitter:
c = temp(index_y,index_x);

%     function z = c(ind_x,ind_y)                
%         if (abs(center_y - ind_y) < 3) && (abs(center_x - ind_x) < 3)
%             z = 98.8e-2;
%         else
%             z = 20e-2;
%         end
%     end

% Anfangsverteilung mit Randwerten = 0
%u(2:mx-1,2:nx-1) = 3 * exp(-0.25*(inner_x.^2+inner_y.^2));

%a = -2;
%b = 5;
%u(2:mx-1,2:nx-1) = a + (b-a).*rand(mx-2,nx-2);

u(2:mx-1,2:nx-1) = random('norm',2,1,[mx-2,nx-2]);
u(floor(center_y)-5:ceil(center_y)+5,floor(center_x)-5:ceil(center_x)+5) = 5*ones(12,12);


% Simulations-figure
figure('Position', [0, 0, 1000, 600]);
handle = surf(x,y,u);
axis manual;
%axis([x_min-border, x_max+border, y_min-border, y_max+border, -1, 1]);
axis([x_min, x_max, y_min, y_max, -1, 5]);



%% Simulation
tic

for i = 1 : 1/dt
%while true 

%      for i = 2 : mx-1
%        for j = 2 : nx-1
%            U(i,j) = dt/h^2 * (u(i-1,j)+u(i+1,j)+u(i,j-1)+u(i,j+1)-4*u(i,j)) + u(i,j);
%        end
%      end
     
     for m = 2 : mx-1
       for n = 2 : nx-1
           U(m,n) = dt/h^2 * (c(m-1,n)*u(m-1,n)+c(m+1,n)*u(m+1,n)+c(m,n-1)*u(m,n-1)+c(m,n+1)*u(m,n+1)-4*c(m,n)*u(m,n)) + u(m,n);
       end
     end
     
    %U = 2*dt*4*c.^2.*del2(u,h,h)+u;
    
    %U(floor(center_y):ceil(center_y),floor(center_x):ceil(center_x)) = U(floor(center_y):ceil(center_y),floor(center_x):ceil(center_x)) + 0.1;
    
    u = U;
    
    % Neumannrand:
    u(:,1) = u(:,2);
    u(1,:) = u(2,:);
    u(mx,:) = u(mx-1,:);
    u(:,nx) = u(:,nx-1);

    
%     % Dirichletrand:
%     u(:,1) = 0;
%     u(1,:) = 0;
%     u(mx,:) = 0;
%     u(:,nx) = 0;
%     
    set(handle, 'ZData', u);
    drawnow
end

toc

disp('Ende')

end