function [theta_pg,phi_pg] = relax_PG(y,Nr,Nt,theta_ini,phi_ini)
eta=1e-7;
tr         =  (0:Nr-1)';
tt         =  (0:Nt-1)';
build_ar   =  @(fr) exp(1i*tr*fr);
build_at   =  @(ft) exp(1i*tt*ft);
theta_pg = theta_ini;
phi_pg = phi_ini;
for i=1:20
    
bb1 = y* conj(build_at(phi_pg));
bb2 = y.'* conj(build_ar(theta_pg));

bb1_re = real(bb1);
bb1_im = imag(bb1);
bb2_re = real(bb2);
bb2_im = imag(bb2);

cos_vec1 = cos(tr.'*theta_pg).';
sin_vec1 = sin(tr.'*theta_pg).';
cos_vec2 = cos(tt.'*phi_pg).';
sin_vec2 = sin(tt.'*phi_pg).';

zz_r1 = cos_vec1.* bb1_re;
zz_r2 = sin_vec1.* bb1_im;
zz_r3 = cos_vec1.* bb1_im;
zz_r4 = sin_vec1.* bb1_re;
grad1_theta = 2*sum(zz_r1+zz_r2)*sum( tr .* (-zz_r4 +zz_r3));
grad2_theta = 2*sum(zz_r3-zz_r4)*sum( tr .* (-zz_r2 -zz_r1));
grad_theta =- grad1_theta -grad2_theta;

zz_t1 = cos_vec2.* bb2_re;
zz_t2 = sin_vec2.* bb2_im;
zz_t3 = cos_vec2.* bb2_im;
zz_t4 = sin_vec2.* bb2_re;
grad1_phi = 2*sum(zz_t1+zz_t2)*sum(tt.* (-zz_t4 +zz_t3));
grad2_phi = 2*sum(zz_t3-zz_t4)*sum(tt.* (-zz_t2 -zz_t1));
grad_phi = -grad1_phi -grad2_phi;
obj_buff =- abs(build_ar(theta_pg)'*bb1)^2;
while 1
        y_theta = theta_pg - eta * grad_theta;
        y_phi = phi_pg - eta*grad_phi;
        At_iter = build_at(y_phi);
        Ar_iter = build_ar(y_theta);
        obj = - norm(Ar_iter'*y*conj(At_iter))^2;
        obj_up = obj_buff  - eta/2*(abs(grad_theta)^2+abs(grad_phi)^2);
        if obj<obj_up 
            
            break
        else
            eta = eta/2;
        end
        if eta<1e-9
            break
        end
end
    theta_pg = y_theta;
    phi_pg = y_phi;
end