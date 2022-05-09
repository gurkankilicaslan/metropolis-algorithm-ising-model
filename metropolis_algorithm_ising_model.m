clear
N = 200; Q = 100000; K = 100; J = -1; kbT = 1.5;
Tmin = 1.5; Tmax = 3; dT = 0.001;
a = ceil(rand(N,N)*2)*2 - 3;
B = 0; % 0, 0.1, 0.5, 1.5
for q = 1 : Q
    kbT = kbT + dT;
    if (kbT < Tmin) || (kbT > Tmax)
        dT = -dT;
    end
    kbT;
        
                  
    r0 = ceil(rand(N^2,2)*N); rn = mod(r0 - 2,N) + 1; rp = mod(r0,N) + 1;
    r = rand(N^2,1);
    
    total_energy = 0;
    spin_total = 0;
    for n = 1 : N^2
        m = (a(rn(n,1),r0(n,2)) + a(rp(n,1),r0(n,2)) + ...
             a(r0(n,1),rn(n,2)) + a(r0(n,1),rp(n,2)))... % sum of outer shell
            *a(r0(n,1),r0(n,2)); % multiply with center divide by 2 add 3
        
        dE= 2*B-2*m;
        if r(n) < exp(-dE/kbT)
            a(r0(n,1),r0(n,2)) = -a(r0(n,1),r0(n,2)); % flip the spin
        end
        total_energy = total_energy + -J*((a(rn(n,1),r0(n,2)) + a(rp(n,1),r0(n,2)) + ...
             a(r0(n,1),rn(n,2)) + a(r0(n,1),rp(n,2)))... % sum of outer shell
            *a(r0(n,1),r0(n,2))) - B*a(r0(n,1),r0(n,2));
        spin_total = spin_total + a(r0(n,1),r0(n,2));
    end
    magnetization = spin_total/N^2;
    disp("kbT total energy and magnetization")
    kbT
    total_energy  
    magnetization 
    
    imagesc(a);
    title(['Monte Carlo Run for Current Temperature']);
    axis equal off;
    drawnow;
end