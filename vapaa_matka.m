function [ la ] = vapaa_matka( temp )

    su = 110.4;
	temp_r = 293.15;
	l_r = 66.5e-9;
	
	eka = temp./temp_r;
	tokayla = 1.0+su/temp_r;
	tokaala = 1.0+su./temp;
	toka = tokayla./tokaala;
	
	la= l_r.*eka.*toka;
    
end

