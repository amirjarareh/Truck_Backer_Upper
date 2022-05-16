%*******************************************
%           AmirHossein Jarareh
%                400616004
%*******************************************
close all;
clear all
clc;

%-------config-------
x_init = 7;
y_init = 10;
phy_init = -90; %deg
n = 300;


%------select mode-----
%1:Generate data
%2:Create table and rules
%3:Prediction Truck
mode = 3; 

if mode == 1
    %Generate data
    [x_g,y_g,ph_g,theta_g,t_g]=generate_raw_data(n,x_init,y_init,phy_init);
    %sort x_g and ph_g and theta_g to X,Phi,Theta
    X = x_g;
    Y = y_g;
    Theta = theta_g;
    
elseif mode == 2
    %Create table and rules
    [ph_num , x_num , theta_num] = lookup_table(Phi,X,Theta);
    ph_num = ph_num';
    x_num = x_num';
    theta_num = theta_num';
    Table = [x_num ph_num theta_num];
    Table_f = ids_filter(Table);
    
elseif mode == 3
    %Prediction Truck
    Table_f = modify_Table;
    [x_out,y_out,phy_out,theta_out,t,ex,eph]=close_loop(n,x_init,y_init,phy_init,Table_f);
    
    figure
    p1 = plot(x_out,y_out,'-.');
    p1(1).LineWidth = 2;
    title('Truck');
    xlabel('x');
    ylabel('y');
    legend('XY');
    grid on
    
    figure
    p1 = plot(t,x_out);
    p1(1).LineWidth = 2;
    title('t and x');
    xlabel('t');
    ylabel('x');
    legend('t,x(t)');
    grid on
    
    figure
    p1 = plot(t,y_out);
    p1(1).LineWidth = 2;
    title('t and y');
    xlabel('t');
    ylabel('y');
    legend('t,y(t)');
    grid on
    
    figure
    p1 = plot(t,phy_out);
    p1(1).LineWidth = 2;
    title('t and phy');
    xlabel('t');
    ylabel('phy');
    legend('t,phy(t)');
    grid on
    
    figure
    t = (1:length(theta_out));
    p1 = plot(t,theta_out);
    p1(1).LineWidth = 2;
    title('t and theta');
    xlabel('t');
    ylabel('theta');
    legend('t,theta(t)');
    grid on

end
% ---------------------All Functions----------------

%generate data
function [x_out,y_out,phy_out,theta_out,t]=generate_raw_data(num,x0,y0,phy0)
b = 4;
x(1) = x0;
y(1) = y0;
phy(1) = phy0 * pi /180;
theta(1) = 0 ;
t = (1:num+1);
i = 1;
while true
        theta_input = input("theta ra vared konid = ");
    if theta_input == 1000
        break
    end
    theta_input = theta_input*pi/180 ;
    theta(i+1) = theta_input;
   
    x(i+1) = x(i) + cos(phy(i)+theta(i+1)) + sin(theta(i+1))*sin(phy(i));
    y(i+1) = y(i) + sin(phy(i)+theta(i+1)) - sin(theta(i+1))*cos(phy(i));
    phy(i+1) = phy(i)-asin(2*sin(theta(i+1))/b);
    plot(x,y)
    grid on
   i=i+1;
end
x_out = x;
y_out = y;
phy_out = phy;
theta_out = theta;
end

%triangle function
function [y] = triangle_func(i,x,a,b,n)
h = (b-a)/n;
if i > n+1
    disp("i is bigger than usual");
elseif i<0 
    disp("i must positive");
else
    if i == 1
        y = triangle(x,a,a,a+h);
    elseif i == n+1
        y = triangle(x,b-h,b,b);
    else
        b1 = a + (i-1)*h;
        a1 = b1 - h;
        c1 = b1 + h;
        y = triangle(x,a1,b1,c1);
    end
end
end
function [y] = triangle(x,a,b,c)
if a == b
    if x <=c && x>=a
        m1 = 1/(c-a);
        y = -1*m1*(x-c);
    else
        y = 0;
    end
elseif b == c
    if (x<=b) && (x>=a)
        m1 = 1/(c-a);
        y = m1*(x-a);
    else
        y = 0;
    end
else
    m = 1/(b-a);
    if (x>=a) && (x<=b)
        y = m*(x-a);
    elseif (x>b) && (x<=c)
        y = -1*m*(x-c);
    else
        y = 0;
    end
end
end

%build lookup table
function [eph_num , ex_num , theta1_num] = lookup_table(eph,ex,theta1)
    %eph_num % % 
    a2 = -90*pi/180;
    b2 = 270*pi/180;
    N2 = 7;
    k2 = 0;
    %ex_num 
    a1 = 0;
    b1 = 20;
    N1 = 5;
    k1 = 0;
    %theta1_num
    a3 = -40*pi/180;
    b3 = 40*pi/180;
    N3 = 7;
    k3 = 0;
    for i = (1:length(ex))
        k1 = 0;
    for j = (1:N1)
        y1 = triangle_func(j,ex(i),a1,b1,N1);
        if y1 ~= 0    
            if y1>=k1
                k1 = y1;
               ex_num(i) = j;
            end
        end
    end
    end

    for i = (1:length(eph))
        k2 = 0;
        k3 = 0;
        for j = (1:N2)
            y2 = triangle_func(j,eph(i),a2,b2,N2);
            if y2 ~= 0    
                if y2>=k2
                    k2 = y1;
                    eph_num(i) = j;
                end
            end
            y3 = triangle_func(j,theta1(i),a3,b3,N3);
            if y3 ~= 0    
                if y3>=k3
                    k3 = y3;
                    theta1_num(i) = j;
                end
            end    
        end
    end
end
function [id_f] = ids_filter(ids_raw)
    i = 0;
    [n,m] = size(ids_raw);
    while i < n
        i = i + 1;
            n_sel = 0;
            selected = 0;
            if i ~= n
                for j = (i+1:n)
                    if ids_raw(i,(1:m-1)) == ids_raw(j,(1:m-1))
                        n_sel = n_sel + 1;
                        selected(n_sel) = j;
                    end
                end  
                if selected ~= 0
                    ids_raw(selected,:) = []; 
                end 
                [n,m] = size(ids_raw);
            end
    end
    id_f = ids_raw;
end


%model system
function [x_out,y_out,phy_out,theta_out,t,ex_out,eph_out]=close_loop(num,x0,y0,phy0,Table_f)
b = 4;
x(1) = x0;
y(1) = y0;
ex(1) = 0;
eph(1) = 0;
phy(1) = phy0*pi/180;
theta(1) = 0 ;
t = (1:num+1);

x_d = 10;
phy_d = 90;

i = 1;
while true
    ex(i) = x_d - x(i);
    eph(i) = phy_d*pi/180 - phy(i) ;
    
    theta(i) = prediction(conv_ex(ex(i)),conv_eph(eph(i)),Table_f);
    x(i+1) = x(i) + cos(phy(i)+theta(i)) + sin(theta(i))*sin(phy(i));
    y(i+1) = y(i) + sin(phy(i)+theta(i)) - sin(theta(i))*cos(phy(i));
    phy(i+1) = phy(i)-asin(2*sin(theta(i))/b);
    
    i=i+1;
    clc
        disp(['predection : ',num2str(round(i/num*100)),'%'])
    if i == num
        break
    end
end
[x_out,y_out,phy_out,theta_out] = map2out(x,y,phy,theta);
ex_out =x_d -x_out;
eph_out =phy_d-phy_out;
t=[1:length(x_out)];
end

%prediction
function [theta_out] = prediction(ex,eph,Table_f)
sum1 = 0 ;
sum2 = 0;

[n,m]=size(Table_f);
a1 = 0;
b1 = 20;
N1 = 5;

a2 = -90*pi/180;
b2 = 270*pi/180;
N2 = 7;

a3 = -40;
b3 = 40;
N3 = 7;
percent = 0;
for i=(1:n)
    row = Table_f(i,:);
    sum1 = sum1 + triangle_func(row(1),ex,a1,b1,N1)* triangle_func(row(2),eph,a2,b2,N2)*(a3 + (row(3)-1)*((b3-a3)/N3))*pi/180;
    sum2 = sum2 + triangle_func(row(1),ex,a1,b1,N1)* triangle_func(row(2),eph,a2,b2,N2);
end
theta_out = sum1/sum2;

end


%data analyze
function [x]=conv_ex(ex)
    x=10-ex;
end
function [phy]=conv_eph(eph)
    phy=pi/2-eph;
end
function [x_out,y_out,phy_out,theta_out] = map2out(x,y,phy,theta)
    x_out = x +0.66;
    y_out = y;
    phy_out = phy*pi/180;
    theta_out = theta*pi/180;
end
function [Table_out] = modify_Table()
    x = get_x();
    phy = get_phy();
    theta = get_theta();
    Table_out = [x' phy' theta'];
end
function [x_out] = get_x()
x_out = [1 1 1 1 1 2 2 2 2 2 2 3 3 3 3 3 4 4 4 4 4 4 5 5 5 5 5 5 5 4 3 2];
end
function [phy_out] = get_phy()
phy_out = [1 2 3 4 5 1 2 3 4 5 6 2 3 4 5 6 2 3 4 5 6 7 3 4 5 6 7 2 1 1 1 7];
end
function [theta_out] = get_theta()
theta_out = [2 2 3 5 7 1 1 2 4 2 3 1 2 4 6 7 5 1 2 5 5 6 1 2 4 3 5 3 2 2 2 1];
end
