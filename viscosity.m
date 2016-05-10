function out = viscosity( temp )
mu0=1.716e-5;
T0=273.11;
S=110.56;

out=mu0*((T0+S)/(temp+S))*(temp/T0)^1.5;

end

