function n = nSiO2w(wl)
    ww = wl/1000;
    a = 0.6961663;
    b = 0.0684043;
    c = 0.4079426;
    d = 0.1162414;
    e = 0.8974794;
    f = 9.896161;
    n = sqrt(1 + a*ww^2/(ww^2-b^2) + c*ww^2/(ww^2-d^2) + e*ww^2/(ww^2-f^2)); 
end