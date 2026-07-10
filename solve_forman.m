% -------------------------------------------------------------------------
% setup system and solve problem
% -------------------------------------------------------------------------
disp('>> setup system and solve problem')
z = fcmplx(1).bc(:,3); 
x = fcmplx(1).bc(:,1); 
y = fcmplx(1).bc(:,2);

A = fCBS_1 * k * fBO_1; 

int_indx = find(z ~= maxx3 & z ~= minx3);

t_old = zeros(fcmplx(1).num(1).val,1); 
t_old(z == maxx3) = 1;
rhs = -A * (z == maxx3);
t_new = A(int_indx, int_indx) \ rhs(int_indx);
t_old(int_indx) = t_new;

flux = abs(fs2 * k * fBO_1 * t_old);

small = 10^-8;

zsf_v = zeros(size(fcmplx(1).bc(:,1)));  
zsf_v(fcmplx(1).bc(:,1) < minx3 + small) = 1;  

zsf_e = abs(fBO_1) * zsf_v;  
zsf_e(zsf_e > 1) = 0;

zflux = zsf_e .* flux;

disp(((sum(zflux) * 10^(-18) * 0.001) / ((maxx3 - minx3) * 10^(-6))))
Kx = ((sum(zflux) * 10^(-18) * 0.001) / ((maxx3 - minx3) * 10^(-6)));

