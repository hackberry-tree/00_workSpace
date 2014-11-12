#!/usr/bin/perl 
# ver Jan2012
use warnings;

sub main {
	my @ld = glob ('*');
	#print "@ld\n";
	system("rm -f sum.out");
	system("./pos-sigma.sh");
	open (SUM, ">> sum.out"); # && print "output : sum.out\n";
        print SUM "# PBE   2a 4f 8i1 8i2 8j :    energy (eV) | cell-parameters    a  b  c (A)     Vol (A³) |   x_4f       x_8i1      y_8i1      x_8i2      y_8i2      x_8j      z_8j  | Mag Mom\n";
	foreach $d (@ld)
	{
		if ((-d $d) && ($d ne "steps"))
		{
			print "FILE $d : \n";
			#print SUM "$d :";
			chdir "$d";
			&dirCheck ();
			chdir "..";	
		}
	}
	close SUM;
}

sub dirCheck {
	my $fe = 0;
	my $vol = 0;
        my $a = 0;
        my $b = 0;
        my $c = 0;
	#print "@ls\n";	
			open (OUT,"OUTCAR");      #here 1
			my $i = 0;
			while (<OUT>)
			{		
                        $i++ if ($i > 0);
                        $i = 0 if ($i == 6);
                        if ($_ =~ m/free  energy   TOTEN/)
                        {
                                my @l1 = split ();
                                $fe = $l1[4];
                        }
                        if ($_ =~ m/volume of cell :/)
                        {
                                my @l2 = split ();
                                $vol = $l2[4];
                                $i = 1;
                        }
                        if ($i == 3)
                        {
                                my @l3 = split ();
                                $a = $l3[0];
                        }
                        if ($i == 4)
                        {
                                my @l4 = split ();
                                $b = $l4[1];
                        }
                        if ($i == 5)
                        {
                                my @l5 = split ();
                                $c = $l5[2];
                        }
                }
                close OUT;
		open (INI,"POSCAR.ini");
                while (<INI>)
                {  
                        if ($_ =~ m/sigma/)
                        {
                                my @ligne = split ();
                                $SN = $ligne[0];
                                $SA = $ligne[1];
                                $SB = $ligne[2];
                                $SC = $ligne[3];
                                $SD = $ligne[4];
                                $SE = $ligne[5];
                        }
                }
                close INI;

               open (OSZ,"OSZICAR");
               $MM = 0;
               while (<OSZ>)
                {
                        if ($_ =~ m/mag/)
                        {
                                my @ligne = split ();
                                $MM = $ligne[9]/30;
                        }
                }
                close OSZ;

	        &anaPOS ();
                my $XB= sprintf "%.8f", $XB; 
                my $XC= sprintf "%.8f", $XC; 
                my $YC= sprintf "%.8f", $YC; 
                my $XD= sprintf "%.8f", $XD; 
                my $YD= sprintf "%.8f", $YD; 
                my $XE= sprintf "%.8f", $XE; 
                my $ZE= sprintf "%.8f", $ZE; 
                print SUM "$SN $SA $SB $SC $SD $SE : $fe eV  $a  $b  $c  $vol $XB $XC $YC $XD $YD $XE $ZE $MM\n";
}



sub anaPOS {
         $XA=0; $YA=0; $ZA=0 ;
         $XB=0; $YB=0; $ZB=0 ;
         $XC=0; $YC=0; $ZC=0 ;
         $XD=0; $YD=0; $ZD=0 ;
         $XE=0; $YE=0; $ZE=0 ;
                $I=0;
		open (POS,"../pos.out");
                while (<POS>)
                {  
                        if ($_ =~ m/position/)
                        {
                                my @ligne = split ();
                                if ($d==$ligne[3])
				{	
                                $I=0;
                                #$FN = $ligne[3];
                                }
                        }
                        if ($I==1)
			{
                                my @ligne = split ();
                                $XA = $ligne[0] ; $YA = $ligne[1] ; $ZA =$ligne[2]; 
			}
                        if ($I==2)
			{
                                my @ligne = split ();
                                $XB = $ligne[0] ; $YB = $ligne[1] ; $ZB =$ligne[2]; 
			}
                        if ($I==3)
			{
                                my @ligne = split ();
                                $XC= $ligne[0] ; $YC = $ligne[1] ; $ZC =$ligne[2]; 
			}
                        if ($I==4)
			{
                                my @ligne = split ();
                                $XD = $ligne[0] ; $YD = $ligne[1] ; $ZD =$ligne[2]; 
			}
                        if ($I==5)
			{
                                my @ligne = split ();
                                $XE = $ligne[0] ; $YE = $ligne[1] ; $ZE =$ligne[2]; 
			}
                       $I++;
                }
                close POS;
}


&main ();
