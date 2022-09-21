%% Lambert's Solver Short Way

r_1= [-6991 -5118 7078]; %First position, in km 
r_2= [11326 -6317 -6677]; %Second Position, in km

delta_t= 5200; %Time between first and second position, in s

R1= norm(r_1); %Distance of the first position, in km 
R2= norm(r_2); %Distance of the second position, in km

R12= R1+R2; %Sum of the distances, in km 
c= norm(r_1-r_2); %Length of chord connecting the distances i.e the absolute value of the difference in length, in km

a_min= (R12 + c)/4; %Minimum semi-major axis, in km

x0 = [a_min 13000]; % initial interval
options = optimset('Display','iter'); % show iterations
[x, fval, exitflag, output] = fzero(@f,x0,options)

function y= f(a_guess)
mu= 3.986*10^5; %Gravitational parameter of Earth, in km^3/s^2

r_1= [-6991 -5118 7078]; %First position, in km 
r_2= [11326 -6317 -6677]; %Second Position, in km

delta_t= 5200; %Time between first and second position, in s

R1= norm(r_1); %Distance of the first position, in km 
R2= norm(r_2); %Distance of the second position, in km

R12= R1+R2; %Sum of the distances, in km 
c= norm(r_1-r_2); %Length of chord connecting the distances i.e the absolute value of the difference in length, in km

a_min= (R12 + c)/4; %Minimum semi-major axis, in km

beta_c= 2*asin(sqrt((R12-c) / (R12+c))); %Critical beta angle, in radians 

if (beta_c/2 <= pi/2)
    b_c= beta_c; %If beta_c/2 in the first quadrant then b_c short = b_c
else 
    b_c= -beta_c; %If beta_c/2 in the fourth quadrant then b_c short = -b_c
end   

n_c= sqrt(mu/a_min^3); %Critical mean anomaly in s^-1
t_c= (1/n_c) * (pi - (b_c - sin(b_c))); %Critical time, in s

alpha_ast= 2*asin((R12 + c) / (4*a_guess)); %Alpha astericks, in radians

if (delta_t < t_c)
    alpha= alpha_ast; %Alpha, in radians
elseif (delta_t > t_c)
    alpha= 2*pi - alpha_ast; %Alpha, in radians
else
    alpha= pi; %Alpha, in radians
end

beta_ast= 2*asin((R12 - c) / (4*a_guess)); %Beta astericks, in radians

beta= beta_ast; %Since going the short way beta= beta_ast, in radians

t= sqrt(a_guess^3/mu) * ( (alpha - sin(alpha)) - (beta - sin(beta)) );

y= delta_t  - t;
end

