var Pie Y R;

varobs Pie Y R;

svar_identification;
exclusion lag 0;
equation 1, R;
exclusion lag 1;
equation 1, Y, R;
equation 2, Pie, R;
equation 3, Pie, Y;
exclusion lag 2;
equation 1, Y, R;
equation 2, Pie, R;
equation 3, Pie, Y, R;
restriction equation 1, coeff(R,0) - coeff(Pie,0) = 0;
end;

svar_global_identification_check;

