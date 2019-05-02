function [tau,r_n] = synch_temp(rl,Fse)

 Ns = length(rl);
 mu =10;
 tau  = 0;
 rle = [zeros(Fse,1); rl; zeros(Fse,1)];
 int_tau = round(tau);
 te = Fse + (1:Fse:(1+(Ns-1)*Fse)) + int_tau;
 frac_tau = tau - int_tau;
 r_n = rle(te)*(1-frac_tau) + rle(te + 1)*frac_tau;
 r_nd = 0.707*(sign(real(r_n)) + 1i * sign(imag(r_n)));
 err = r_n - r_nd;
 drl = 0.5*(1 - frac_tau)*(rle(te+1)-rle(te-1)) + 0.5*frac_tau*(rle(te+2)-rle(te));
 tau = tau - mu * real(err'*drl/Ns);
end