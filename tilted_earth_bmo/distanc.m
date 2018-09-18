function d=distanc(k,l)
% fprintf('in fun\n');
% disp(k);
% fprintf('in fun\n');
[p1,t1,r1] = cart2sph(k(1),k(2),k(3));
[p2,t2,r2] = cart2sph(l(1),l(2),l(3));
hsi = ((1-cos(t1-t2))./2)+((cos(t1))*(cos(t2))*((1-cos(p1-p2))./2));
d = 2*r1*asin(sqrt(hsi));
end
