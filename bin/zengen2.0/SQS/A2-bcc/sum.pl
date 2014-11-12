#!/usr/bin/perl 
use warnings;

sub main {
	my @ld = glob ('*');
	system("rm -f sum.out");
	open (SUM, ">> sum.out"); # && print "output : sum.out\n";
        print SUM "# x_B    energy      | cell-param    a   c   Vol  | Mag Mom \n";

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
			open (OUT,"OUTCAR");
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

               open (OSZ,"OSZICAR");
               $MM = 0;
               while (<OSZ>)
                {
                        if ($_ =~ m/mag/)
                        {
                                my @ligne = split ();
                                $MM = $ligne[9]/16;    # FOR 16 ATOMS !!!
                        }
                }
                close OSZ;

                printf SUM "%4s : %10.4f eV | %8.4f %8.4f %8.4f |  %8.4f muB/at. \n", $d, $fe, $a, $c, $vol, $MM ;
}

&main ();
