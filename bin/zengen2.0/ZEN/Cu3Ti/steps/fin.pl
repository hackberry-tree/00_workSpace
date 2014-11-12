#!/usr/bin/perl
#;-*- Perl -*-
# code made by R.Souques 2011
use warnings;

sub main {
	my @ld = glob ('*');
	#print "@ld\n";
	system("rm -f brut.out");
	open (SUM, ">> brut.out"); 
	for ($d = 1; $d <= $#ld; $d++)
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
        &concatefiles
}

sub dirCheck {
	my $fe = 0;
	my $vol = 0;
	my @ls = glob ('*');
	#print "@ls\n";	
	foreach $f (@ls)
	{
		if ($f =~ m/OUTCAR/)
		{
			open (OUT,"$f");
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
			}
			print "$f ; ";
			print SUM " $vol $fe\n";
			close OUT;	
		}
	}
	print "\n";	
}

 sub concatefiles {
	system("rm -f results.out");
	open (RES, ">> results.out"); 
        @ARGV = ("conf.out","brut.out");
        while ($ligne = <>)
          {
          print SUM "$ligne\n";
          }
        close RES;
}

&main ();
