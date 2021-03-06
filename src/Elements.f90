Module Elements
Implicit Real(8) (A-H,O-Z)

Character(2) ElName(120)
Real(8)		AMS(120),		& ! Isotope atomic masses
			RA(120)			  ! Covalent radii by Melvin-Hughes


      DATA ELNAME/ &
     'H ','He', &
     'Li','Be','B ','C ','N ','O ','F ','Ne', &
     'Na','Mg','Al','Si','P ','S ','Cl','Ar', &
     'K ','Ca',								  &
     'Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu','Zn', &
               'Ga','Ge','As','Se','Br','Kr',			&
     'Rb','Sr',											&
     'Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd',	&
               'In','Sn','Sb','Te','I ','Xe',			&
     'Cs','Ba',											&
     'La',14*'XY','Hf','Ta','W ','Re','Os','Ir','Pt','Au','Hg', &
               'Tl','Pb','Bi','Po','At','Rn',11*'  ','Xx','Du',21*'  '/

	    DATA (AMS(Iams),Iams=1,54)  / & !     Isotope masses
        1.007825D+00,4.0026D+00,7.01600D+00,9.01218D+00,11.00931D+00, &
        12.0D+00,14.00307D+00,15.99491D+00,18.99840D+00,19.99244D+00, &
        22.9898D+00,23.98504D+00,26.98153D+00,27.97693D+00,			  &
        30.97376D+00,31.97207D+00,34.96885D+00,39.948D+00,			  &
        38.96371D+00,39.96259D+00,44.95592D+00,47.90D+00,50.9440D+00, &
        51.9405D+00,54.9381D+00,55.9349D+00,58.9332D+00,57.9353D+00,  &
        62.9298D+00,63.9291D+00,68.9257D+00,73.9219D+00,74.9216D+00,  &
        79.9165D+00,78.9183D+00,83.9115D+00,						  &
        84.9117D+00,87.9056D+00,89.9054D+00,89.9043D+00,92.9060D+00,  &
        97.9055D+00,97.0D+00,101.9037D+00,102.9048D+00,105.9032D+00,  &
        106.9041D+00,113.9036D+00,114.9041D+00,119.9022D+00,		  &
        120.9038D+00,129.9067D+00,126.9044D+00,131.9042D+00/		  
      DATA (AMS(Iams),Iams=55,106)  /										  &
        132.9054D+00,137.9052D+00,138.9063D+00,139.9054D+00,		  &
        140.9076D+00,141.9077D+00,144.9127D+00,151.9197D+00,		  &
        152.9212D+00,157.9241D+00,158.9253D+00,163.9292D+00,		  &
        164.9303D+00,165.9303D+00,168.9342D+00,173.9389D+00,		  &
        174.9408D+00,179.9465D+00,180.9480D+00,183.9509D+00,		  &
        186.9557D+00,191.9615D+00,192.9629D+00,194.9648D+00,		  &
        196.9665D+00,201.9706D+00,									  &
        204.9744D+00,207.9766D+00,208.9804D+00,208.9824D+00,		  &
        209.9871D+00,222.0176D+00,									  &
        223.0197D+00,226.0254D+00,									  &
        227.0278D+00,232.0381D+00,231.0359D+00,238.0508D+00,		  &
        237.0482D+00,244.0642D+00,243.0614D+00,247.0703D+00,		  &
        247.0703D+00,251.0796D+00,252.0829D+00,257.0751D+00,		  &
        258.0986D+00,259.1009D+00,260.1053D+00,261.1087D+00,		  &
        2*0.0D+00/
      DATA (AMS(Iams),Iams=107,120)  /14*0.d0/


     DATA RA/ &        !RA is a covalent atom radius by Melvin-Hughes. Unknown radii taken as 1.000
     0.3707,0.53,&
     1.520,1.113,0.795,0.771,0.547,0.70006537,0.709,1.60,&
     1.858,1.599,1.432,1.176,0.947,1.02,0.994,1.92,	 &
     2.272,1.974,1.00000,1.44,1.32,1.249,1.366,1.241,1.253,1.246,&
                 1.278,1.333,1.221,1.149,1.248,1.16,1.1415,1.97, &
     2.475,2.151,1.00000,1.59,1.00000,1.363,1.00000,1.34,1.35,1.376,&
                 1.445,1.49,1.626,1.405,1.45,1.35,1.3333,2.18,		&
     2.655,2.174,1.870,1.825,9*1.00000,1.86,3*1.00000,1.61,1.00000,	&
                             1.371,1.00000,1.338,1.00000,1.388,		&
                 1.442,1.503,1.704,1.75,1.548,8*1.00000,1.49,5*1.00000,0.1,1.0,&
     21*1.d0/


End Module

