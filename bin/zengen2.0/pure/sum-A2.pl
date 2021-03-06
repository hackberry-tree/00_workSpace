#!/usr/bin/perl 
use warnings;

sub main {
	my @ld = glob ('*');
	#print "@ld\n";
	system("rm -f sum.out");
	open (SUM, ">> sum.out"); # && print "output : sum.out\n";
        print SUM "## VASPv5 - PBE - 400ev - 21x21x21 - NBANDSok - no SP\n"; 
        print SUM "##          Total-energy | mag=   |  Vol-at  Vol-cub  a-cub \n"; 
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
			open (OUT,"OUTCAR.400");
			my $i = 0;
			while (<OUT>)
			{		
                        $i++ if ($i > 0);
                        $i = 0 if ($i == 6);
                        if ($_ =~ m/free  energy   TOTEN/)
                        {
                                my @l1 = split ();
                                $fe = $l1[4]
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
                $MO=0;
                open (OSZ,"OSZICAR.400");
                while (<OSZ>)
                {
                        if ($_ =~ m/mag/)
                        {
                                my @ligne = split ();
                                $MO = $ligne[9];
                        }
                }
                close OSZ;
                my $acub=$a*2;
                my $vcub=$acub**3;
                my $vat=$vcub/4;
               # print SUM "$d : $fe eV  $vcub angs^3 $acub $mag $MO \n";
                printf SUM "%6s : %12.5f eV | %6.3f | %7.2f %7.2f %6.2f\n", $d, $fe, $MO, $vat, $vcub, $acub;
}

&main ();
