#!/usr/bin/perl 
# ver Jan2012
use warnings;

sub main {
	my @ld = glob ('*');
	#print "@ld\n";
	system("rm -f sum.out");
        system("./pos-chi.sh");
	open (SUM, ">> sum.out"); # && print "output : sum.out\n";
        print SUM "#PBE 2a 8c 24g1 24g2 :        energy | cell-para (Angs)   a   b   c    Volume (A^3) |  X_8c        X_24g1      Z_24g1      X_24g2     Z_24g2  | MAG MOM\n";
	foreach $d (@ld)
	{
		if (-d $d)
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
                        if ($_ =~ m/chi/)
                        {
                                my @ligne = split ();
                                $SN = $ligne[0];
                                $SA = $ligne[1];
                                $SB = $ligne[2];
                                $SC = $ligne[3];
                                $SD = $ligne[4];
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
                                $MM = $ligne[9]/29;
                        }
                }
                close OSZ;

                &anaPOS ();
                my $A= sprintf "%.8f", $XBB;
                my $B= sprintf "%.8f", $XCC;
                my $C= sprintf "%.8f", $ZCC;
                my $D= sprintf "%.8f", $XDD;
                my $E= sprintf "%.8f", $ZDD;
                $a1=-$a*2;$b1=-$b*2;$c1=-$c*2;
                $vol=$vol*2;
                $fe=$fe*2;
                print SUM "$SN $SA $SB $SC $SD  : $fe eV  $a1  $b1  $c1  $vol  $A  $B  $C  $D $E $MM\n";
}

sub anaPOS {
         $XA=0; $YA=0; $ZA=0 ;
         $XB=0; $YB=0; $ZB=0 ;
         $XC=0; $YC=0; $ZC=0 ;
         $XD=0; $YD=0; $ZD=0 ;
         $XAA=0; $YAA=0; $ZAA=0 ;
         $XBB=0; $YBB=0; $ZBB=0 ;
         $XCC=0; $YCC=0; $ZCC=0 ;
         $XDD=0; $YDD=0; $ZDD=0 ;
                $I=0;
                open (POS,"../pos.out");
                while (<POS>)
                {
                        if ($_ =~ m/position/)
                        {
                                my @ligne = split ();
                                if ($d eq $ligne[3])
                                {
                                $I=0;
                                #$FN = $ligne[3];
                                }
                        }
                        if ($I==1)
                        {
                                my @ligne = split ();
                                $XA = $ligne[0] ; $YA = $ligne[1] ; $ZA =$ligne[2];
                                   $XAA=-$XA/2+$YA/2+$ZA/2;
                                   $YAA= $XA/2-$YA/2+$ZA/2;
                                   $ZAA= $XA/2+$YA/2-$ZA/2;
                        }
                        if ($I==2)
                        {
                                my @ligne = split ();
                                $XB = $ligne[0] ; $YB = $ligne[1] ; $ZB =$ligne[2];
                                  $XBB=-$XB/2+$YB/2+$ZB/2;
                                  $YBB= $XB/2-$YB/2+$ZB/2;
                                  $ZBB= $XB/2+$YB/2-$ZB/2;  
                        }
                        if ($I==3)
                        {
                                my @ligne = split ();
                                $XC= $ligne[0] ; $YC = $ligne[1] ; $ZC =$ligne[2];
                                   $XCC=-$XC/2+$YC/2+$ZC/2;
                                   $YCC= $XC/2-$YC/2+$ZC/2;
                                   $ZCC= $XC/2+$YC/2-$ZC/2;
                        }
                        if ($I==4)
                        {
                                my @ligne = split ();
                                $XD = $ligne[0] ; $YD = $ligne[1] ; $ZD =$ligne[2];
                                  $XDD=-$XD/2+$YD/2+$ZD/2;
                                  $YDD= $XD/2-$YD/2+$ZD/2;
                                  $ZDD= $XD/2+$YD/2-$ZD/2;
                        }
                       $I++;
                }
                close POS;
}

&main ();
