function [r_min, den_opt, num_min, s_opt, ts_opt] = find_transfer_params(H, Ts)

r_min = inf;
den_opt = [0,0,0];
num_min = 0;
a1DEL = 10;
a2DEL = 100;
a1LB = 0;
a1UB = 1000;
a2LB = 0;
a2UB = 10000;
s_opt = [];
ts_opt = [];
numLB = min(H);
numUB = max(H);
numDEL = (numUB-numLB)/10;

for LLL=1:10
    if LLL~=1
        a1LB = den_opt(2)-a1DEL;
        a1UB = den_opt(2)+a1DEL;
        a1DEL = a1DEL/10;

        a2LB = den_opt(1)-a2DEL;
        a2UB = den_opt(1)+a2DEL;
        a2DEL = a2DEL/10;
        
        numLB = num_min-numDEL/2;
        numUB = num_min+numDEL/2;
        numDEL = numDEL/10;
    end
    for num = numLB:numDEL:numUB
        for alpha1=a1LB:a1DEL:a1UB
            for alpha2=a2LB:a2DEL:a2UB
                den = [alpha2, alpha1, 1];
                [s, ts] = sisoctf2dstep(num, den, 0, Ts, size(H,1)-1);
                r_new = sum((H-s(1:end)).^2);
                if r_new < r_min
                    r_min = r_new;
                    den_opt = den;
                    num_min = num;
                    s_opt = s;
                    ts_opt = ts;
                end
            end
        end
    end
    disp(LLL);
end


end