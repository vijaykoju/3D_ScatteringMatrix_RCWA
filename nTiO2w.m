function n = nTiO2w(wl)
    ww = wl/1000;
    n = sqrt(5.913 + 0.2441/(ww^2-0.0803));
end