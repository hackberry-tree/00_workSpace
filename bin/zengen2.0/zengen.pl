#!/usr/bin/perl
#;-*- Perl -*-
# code made by R. Souques 2011, updated by J.C. Crivello
# crivello@icmpe.cnrs.fr
#
# 31th May 2013
$ver=2.01;

##########################################
## TO DO LIST
# dft2tdb 
# test INPUT file in ZEN/folders
# case of pure
# fin.pl : ajouter Delta_H
# mixing.pl : trier concentration  

use warnings;
use Cwd;

#######################################################"
#
# where is ZEN folder ?
$zenfolder='/home/enoki/Dropbox/00_scripts/bin/zengen2.0/'; 
$potdir='/home/enoki/vasp/vasp_potential/vasp/potpaw_PBE';

# @El    : table of elements 
# @M     : table of multiplicity
# nM     : number of inequivalent sites
# nC     : number of ordered configurations
# Comb   : table of combinations  

#######################################################"

sub main
{
        system("clear");
        print "############################## ZEN_GEN # version $ver \n";
        # 
        # Ask phase and elements
	&askElPh ();
        # 
        if ((($pha eq 'fcc')||($pha eq 'bcc'))||($pha eq 'hcp'))
        { 
               # 
               # Generate SQS excution files
               $michi="$zenfolder/SQS";
               &sqs ();
               die ("Preparation done\n"); 
               # 
        } 
        $michi="$zenfolder/ZEN";
        if ($El[0] eq 'ALL')
        { 
               # 
               # Cas of pure phase calculated for every element
               &pure ();
               die ("Preparation done\n"); 
        }
        # 
        # Get multiplicity of each Wyckoff site in the considered phase          
	&genM ();
        # 
	$nM = scalar(@M);
	print "#######################################\n";
	print "Studied phase                    = $phase \n";
	print "nb of inequivalent sites         = $nM \n";
	print "table of multiplicity            = @M\n";
	$nEl = scalar(@El);
	print "table of element                   @El\n";
	print "nb of element                    = $nEl \n";
	$nC = $nEl**($nM);
	print "nb of ordered configurations     = $nC \n";
	my @Comb;
	print "#######################################\n";
        print "It is alright? [y/n]\n ";
        $ok = <STDIN>;
        chomp $ok;
        if (($ok eq 'y')||($ok eq 'Y'))
        { print "\n Procedure is starting\n"; sleep(2);}
        else { die ("nothing done.\n"); }
        # 
        # Generate every ordered configuration
	&genC ();
        # 
        # Calculation of concentration and degenerated-concentration (blocposition) of each configuration
	&Xcalc ();
        # 
        # Sort configuration by increasing blocposition, with first: A, B, C... 
	&arrangBP ();  #à virer pour eviter de ralentir
        # 
        # Sort configuration by increasing concentration, with first: uniaries, binaries, ternaries ...
	&arrangX ();
        # 
        # Sort configuration by increasing A, B, C... elements
	&arrangEl ();
        # 
        # Print the final Table of configurations
	&printA ();
        # 
        # Build folders with configuration-coresponding name
	&genF ();
        # 
        # Build a conf.out file which summarizes the Table of configurations
	&genSUM ();
        # 
        # Generate the POSCAR of every configuration
	&genPOS ();
        # 
        # Generate the POTCAR of every configuration
        &genPOT ();
        # 
        # Preparation of execution files
        &preEXE ();
	print ("Generation done.\n");
}
#######################################################"

# Ask phase and elements
sub askElPh
{
	$test_ph = 0;
	$test_el = 0;
	while ($test_ph eq 0) 
	{
	print "#################################################################################\n";
	print "## 1) choice of the studied phase: \n";
	print "##################################################################################\n";
	print " code     # phase   #  prototype   # Space Group    # Pearson  # Nb  # complement \n";
	print "##################################################################################\n";
	print "c (chi)   # A12     #  alpha-Mn    # I-43m    (217) #  cI58    #  4  # chi phase \n";
	print "b (betaMn)# A13     #  beta-Mn     # P4_132   (213) #  cP20    #  2  # beta-Mn \n";
	print "C14       # C14     #  MgZn2       # P6_3/mmc (194) #  hP12    #  3  # Laves AB2 \n";
	print "C15       # C15     #  Cu2Mg       # Fd-3m    (227) #  cF24    #  2  # Laves AB2 \n";
	print "C36       # C36     #  MgNi2       # P6_3/mmc (194) #  hP24    #  5  # Laves AB2 \n";
	print "C36-3s    # C36-3ss #  MgNi2       # P6_3/mmc (194) #  hP24    #  3  # 3 sublattices [A1, A2, B]\n";
	print "s (sigma) # D8b     #  CrFe        # P4_2/mnm (136) #  tP30    #  5  # D8b - sigma phase \n";
	print "d (delta) #         #  Mo7Ni7      # P212121   (19) #  oP56    # 12  # delta-MoNi\n";
	print "delta-sim #         #  Mo7Ni7      # P212121   (19) #  oP56    #  3  # delta simplified [CN12,14,15+]\n";
	print "P         #         #  Cr18Mo42Ni40# Pnma      (62) #  oP56    # 12  # P-CrMoNi \n";
	print "P-sim     #         #  Cr18Mo42Ni40# Pnma      (62) #  oP56    #  3  # P simplified [CN12,14,15+] \n";
	print "Cu3Ti     #         #  Cu3Ti,Ni3Ta # Pmmn      (59) #   oP8    #  3  # beta-Cu3Ti (LT) \n";
	
	print "#################################################################################\n";
	print "# SQS 16  : fcc  bcc  hcp                                  ONLY BINARY SYSTEM     \n";
	print "#################################################################################\n";

		print "Phase? \n";
		chomp ($pha = <STDIN>);
		
        if ($pha eq "s") { $phase = "sigma"; } 
        elsif ($pha eq "c") { $phase = "chi"; }
        elsif ($pha eq "d") { $phase = "delta"; } 
        elsif ($pha eq "ds") { $phase = "delta-sim"; }
        elsif ($pha eq "b") { $phase = "betaMn"; }
        elsif ($pha eq "Ps") { $phase = "P-sim"; }
        elsif ($pha eq "fcc") { $phase = "A1-fcc"; } 
        elsif ($pha eq "bcc") { $phase = "A2-bcc"; }
        elsif ($pha eq "hcp") { $phase = "A3-hcp"; }
        else {$phase = $pha;}
        
        $path_ph = "$zenfolder/ZEN/$phase";
        if ( !$phase ) { next;}
        elsif (-e $path_ph) 
        { 
        	$test_ph = 1;
        	last;
       	}
        else 
        {
        	print "#################################################################################\n";
        	print "!! This phase has no any linked input file !!\n";
        	print "Please add new input files associated to # $pha # in the $zenfolder/ZEN folder,\n";
        	print "or refer to the following phase list:\n";
        	print "#################################################################################\n";
	$test_ph = 0;
        }
	}
	
	print "#################################################################################\n";
	print "## 2) choice of the element(s)?  *** ALL (pure)\n";
	print "#################################################################################\n";
        my $pwd = getcwd(); 
        chdir "$potdir";
        system ("ls -d */");                               # list Potentiels available 
        chdir "$pwd";
	print "#################################################################################\n";
	
	while ($test_el eq 0)
	{
		print "Element(s) separated by a space?  [A B ...] \n";
		chomp (my $tEl = <STDIN>);                         # input of elements without "\n" at the end
		if ( !$tEl ) { next; }
		@El = split (/\s+/,$tEl);                          # build a table of elements from the input
		
		my $n_El = $#El+1;
		$n_test = 0;
		foreach $el (@El)
		{
			$path_el = "$potdir/$el";
	        if (-e $path_el) { ++$n_test; }
	        else 
	        {
	        	print "#################################################################################\n";
                        print "!! Element # $El # has no any linked POTCAR file !!\n";
	        	print "Please add the required POTCAR file in the $zenfolder/POT folder,\n";
	        	last;
	        }
		}
		if ($n_El eq $n_test) 
		{ 
			$test_el = 1;
			last;
		}
		else 
		{
			print "Or, refer to the following element list:\n";
			print "#################################################################################\n";
		        my $pwd = getcwd(); 
		        chdir "$potdir";
		        system ("ls -d */");                               # list Potentiels available 
		        chdir "$pwd";
			print "#################################################################################\n";
		}
	}         
 
}
#######################################################"

# Get multiplicity of each Wyckoff site in the considered phase          
sub genM
{
	@M = ();
	for (my $j = 1 ; $j < 100 ; $j++)
	{
		open (POS, "< $michi/$phase/POS.$j") || last;
		my $n = 0;
		while (<POS>)
		{
			$n++;
		}
		push (@M,"$n");
	}
}
#######################################################"

# Generate every ordered configuration
sub genC
{
	my $j = 0;
        while ($j < $nM) 
	{
                my $i = 0;	
		while ($i < $nC)
		{
			foreach $el (@El)
			{
				for (my $k = 0 ; $k < $nEl**($j) ; $k++)
				{
					$Comb [$i] [$j] = $el;
					$i++;
				}
			}	
		}
		$j++;
	}
}
#######################################################"

# Calculation of concentration and degenerated-concentration (blocposition) of each configuration
sub Xcalc
{
	my $nTot = 0;
	foreach my $t (@M)
	{
		$nTot += $t;
	}
	
	my $p = 0;
	foreach my $el (@El)
	{
		$i = 0;	
		while ($i < $nC)
		{
			my $x = 0;
			for (my $j = 0 ; $j < $nM ; $j++)
			{
				$x += $M[$j] if ($Comb [$i] [$j] eq $el);
			}
			$x = ($x/$nTot);
			$Comb [$i] [$nM+$p] = $x;				#concentration
			$Comb [$i] [$nM+$nEl] += 10+$p  if ($x != 0);		#blocposition
			$i++;
		}
		$p++;
	}
}
#######################################################"

# Sort configuration by increasing blocposition, with first: A, B, C... 
sub arrangBP
{
#	print "1\n";
	my $n = 0;    # nb iterration (just an info)
	my $i = 0;
		while ($i < $nC)
		{
#			print "2\n";
			$k = $i+1;
			while ($k < $nC) 
			{
#				print "3\n";
				if ($k > $i && $Comb[$k][$nM+$nEl] < $Comb[$i][$nM+$nEl])
				{
					($Comb[$i],$Comb[$k]) = ($Comb[$k],$Comb[$i]);	
				}
				$k++;
				$n++;
                                if ( $n % 100 == 0)
			        {
                                        print ".";
                                        if ( $n % 100000 == 0) {print "(1) being sort: please wait";}
			        }
			}
			$i++;
		}
                print "ok \n";
}
#######################################################"

# Sort configuration by increasing concentration, with first: uniaries, binaries, ternaries ...
sub arrangX
{
	my $n = 0;    # nb iterration (just an info)
	my $i = 0;
	while ($i < $nC)
	{
		my $shiftEl = 0;
		for (my $l = 1 ; $l < $nEl ; $l++)
		{
			my $stop = 0;
			my $bEl = $El[$l];
			$shiftEl++; 
			for (my $j = 0 ; $j < $nM ; $j++)
			{
				if ($Comb[$i][$j] eq $bEl)
				{
					$stop = 1;
					last;
				}	
			}
			last if ($stop == 1);
		}
		$k = $i+1;
		while ($k < $nC) 
		{
			if ($Comb[$k][$nM+$nEl] == $Comb[$i][$nM+$nEl] && $Comb[$k][$nM+$shiftEl] < $Comb[$i][$nM+$shiftEl])
			{
				($Comb[$i],$Comb[$k]) = ($Comb[$k],$Comb[$i]);	
			}
			$k++;
			$n++;
                        if ( $n % 100 == 0)
		        {
                                print ".";
                                if ( $n % 100000 == 0) {print "(2) being sort: please wait";}
		        }
		}
		$i++;
	}
        print "ok \n";
}
#######################################################"

# Sort configuration by increasing A, B, C... elements
sub arrangEl
{
	my $nb = 0;    # nb iterration (just an info)
	my $i = 0;
	while ($i < $nC)
	{
		my $shiftEl = 0;
		for (my $l = 1 ; $l < $nEl ; $l++)
		{
			my $stop = 0;
			$bEl = $El[$l];
			$shiftEl++; 
			for (my $j = 0 ; $j < $nM ; $j++)
			{
				if ($Comb[$i][$j] eq $bEl)
				{
					$stop = 1;
					last;
				}	
			}
			last if ($stop == 1);
		}
		my $n = 0;
		for (my $j = 0 ; $j < $nM ; $j++)
		{
			$n += 2*($j+1) if ($Comb[$i][$j] eq $bEl);
		}
		$k = $i+1;
		while ($k < $nC) 
		{
			my $m = 0;
			for (my $j = 0 ; $j < $nM ; $j++)
			{
				$m += 2*($j+1) if ($Comb[$k][$j] eq $bEl);
			$nb++;
                        if ( $nb % 100 == 0)
		        {
                                print ".";
                                if ( $nb % 100000 == 0) {print "(3) being sort: please wait";}
		        }
			}
			if ($m < $n && $Comb[$k][$nM+$nEl] == $Comb[$i][$nM+$nEl] && $Comb[$k][$nM+$shiftEl] <= $Comb[$i][$nM+$shiftEl])
			{
				($Comb[$k],$Comb[$i]) = ($Comb[$i],$Comb[$k]);
				$n = $m;
			}
			$k++;
		}
		$i++;
	}
        print "ok \n";
}
#######################################################"

# Print the final Table of configurations
sub printA
{
	print "#######################################\n";
	my $i = 0;
	while ($i < $nC)
	{
		my $j = 0;
		$n = $i+1;
		if (length ($n) == 1)
		{
			print "$n     ";
		}
		elsif (length ($n) == 2)
		{
			print "$n    ";
		}
		elsif (length ($n) == 3)
		{
			print "$n   ";
		}
		else
		{
			print "$n  ";
		}
				
		while ($j < $nM+$nEl) 
		{
			my $vvv = $Comb [$i] [$j];
			if ($j < $nM)
			{
				if (length ($vvv) == 1)
				{
					print "$vvv ";
				}
				else
				{
					print "$vvv";
				}
			}
			else
			{
				printf '%.3f', $vvv;
			}
			print "  ";
			$j++;
		}
		print "\n";
		$i++;
	}
	print "#######################################\n";
	print " \n";
	print " \n";
}
#######################################################"

# Build folders with configuration-coresponding name
sub genF
{
	print (" writing folders, please wait\n");
	my @name = ();
	push (@name, "$phase");
	foreach my $el (@El)
	{
		push (@name, "$el");
	}
	$file = join ("-", @name);
	system ("mkdir -p $file");
	for (my $i = 1 ; $i <= $nC ; $i++)
	{
		system ("mkdir -p $file/$i");
	}
}
#######################################################"

# Build a conf.out file which summarizes the Table of configurations
sub genSUM
{
	print (" writing files (3), please wait\n");
	open (SUM, "> $file/conf.out");
	my $i = 0;
	while ($i < $nC)
	{
		my $j = 0;
		$n = $i+1;
		if (length ($n) == 1)
		{
			print SUM "$n     ";
		}
		elsif (length ($n) == 2)
		{
			print SUM "$n    ";
		}
		elsif (length ($n) == 3)
		{
			print SUM "$n   ";
		}
		else
		{
			print SUM "$n  ";
		}
				
		while ($j < $nM+$nEl) 
		{
			my $vvv = $Comb [$i] [$j];
			if ($j < $nM)
			{
				if (length ($vvv) == 1)
				{
					print SUM "$vvv ";
				}
				else
				{
					print SUM "$vvv";
				}
			}
			else
			{
				printf SUM '%.3f', $vvv;
			}
			print SUM "  ";
			$j++;
		}
		print SUM "\n";
		$i++;
	}
	close SUM;
}
#######################################################"

# Generate the POTCAR of every configuration
sub genPOT
{
	print (" writing files (1), please wait\n");
	foreach my $el (@El)
        {
		  system ("sh $zenfolder/genpot.sh $el");
        }
        open (SUM, "< $file/conf.out");
	while (<SUM>)
	{
		  my @tabelem = ();
		  my @tabcomp = ();
                  my @line = split ();
                        $folder = $line[0];
                        open (POT, "> $file/$folder/POTCAR");
                        for (my $i = 0 ; $i < $nEl ; $i++)
                        {
                             #$tabelem[$i] = $line[1+$i];
                             $tabcomp[$i] = $line[1+$nM+$i];
		             if ($tabcomp[$i] ne "0.000")
                             {
                                   open (POTEL, "< POTCAR.$El[$i]") || die ("No POTCAR file associated to the $El[$i] element!");
                                   while (<POTEL>)
                                   {
                                        print POT; 
                                   }
                                   close POTEL;
                             }
                        }
                        close POT; 
	}
	close SUM;
        system ("rm POTCAR.*");
}
#######################################################"

# Generate the POSCAR of every configuration
sub genPOS
{
	print (" writing files (2), please wait\n");
	for (my $i = 0 ; $i < $nC ; $i++)
	{
		$f = $i+1;
		system ("cp $michi/$phase/POS.0 $file/$f/temPOS");
		my @nnl = ();
		foreach my $l (@El)
		{
			my $nl = 0;
			for (my $j = 1 ; $j <= $nM ; $j++)
			{
				if ($Comb[$i][$j-1] eq $l)
				{
					open (POS, "< $michi/$phase/POS.$j");
					open (tPOS, ">> $file/$f/temPOS");
					while (<POS>)
					{
						$nl++;
						if (/el/)
						{
							$_ =~ s/el/$l/;
							print tPOS "$_";
						}
						else
						{
							print tPOS "$_";
						}
								
					}
					close (POS);
					close tPOS;
				}
			}
			push (@nnl,"$nl");
		}
		open (tPOS, "< $file/$f/temPOS");
		open (POS, "> $file/$f/POSCAR.ini");
		my $k = 0;
		while (<tPOS>)
		{
			$k++;
			if ($k == 1)
			{
				@cEl = ();
				push (@cEl,"$phase-$f-");
				for (my $z = 0 ; $z < $nM ; $z++)
				{
					my $tEl = $Comb[$i][$z];
					push (@cEl,"$tEl");
				}
				print POS "@cEl\n";
			}
			elsif ($k == 6)
			{
				foreach $toto (@nnl) 
				{
					print POS "  $toto" unless ($toto == 0);
				}
				print POS "\n";
			}
			unless ($k == 1 || $k == 6)
			{		
				print POS "$_";
			}
			
		}
		close (tPOS);
		close (POS);
		system (" rm $file/$f/temPOS");			
	}
}
#######################################################"

# Preparation of execution files
sub preEXE
{
        system ("mkdir -p $file/steps");
        system ("cp $michi/$phase/steps/* $file/steps/.");
        system ("cp $michi/$phase/POSCAR.ini $file/.");
        system ("mv $file/steps/*.pl $file/.");
        system ("mv $file/steps/*.sh $file/.");
}
#######################################################"

# Cas of pure phase calculated for every element
sub pure
{
        my $folder = "PURE-$phase" ;
        system ("mkdir -p $folder");
        system ("mkdir -p $folder/steps");
        system ("cp $michi/$phase/steps/* $folder/steps/.");
        system ("mv $folder/steps/*.pl $folder/.");
        system ("mv $folder/steps/*.sh $folder/.");

        chdir "$folder";
        system ("cp -r $michi/PURE/* .");
        my @ld = glob ('*');
        foreach $d (@ld)
        {
                if (-d $d)
                {
                        print "$d : \n";
                        system ("cp $michi/$phase/POSCAR.ini $d/.");
                }
        }
        chdir "..";
        print "################################################# \n";
        print "# Remplacer: \n";
        print "################################################# \n";
        print "for STEP in `seq \$SMIN \$SMAX`; do \n";
        print "for DIR  in `seq \$CMIN \$CMAX`; do \n";
        print ".......... \n";
        print "done \n";
        print "done \n";
        print "################################################# \n";
        print "# Par: \n";
        print "################################################# \n";
        print "for STEP in `seq \$SMIN \$SMAX`; do \n";
        print "for DIR  in *; do \n";
        print "  if [ -d \$DIR ]; then \n";
        print "  .......... \n";
        print "  fi         \n";
        print "done \n";
        print "done \n";
        print "################################################# \n";

}
#######################################################"

# Generate SQS excution files
sub sqs
{
        $nEl = 2;      # limitation to binary system
        my @name = ();
        push (@name, "$phase");
        for (my $i = 0 ; $i < $nEl ; $i++)
        {
                  system ("sh $zenfolder/genpot.sh $El[$i]");
                  push (@name, "$El[$i]");
        }
        $file = join ("-", @name);
        system ("mkdir -p $file");

        system ("cp -r $michi/$phase/* $file/.");

        open (POT, "> $file/POTCAR");
        for (my $i = 0 ; $i < $nEl ; $i++)
        {
                  open (POTEL, "< POTCAR.$El[$i]");
                  while (<POTEL>)
                  {
                            print POT;
                  }
                  close POTEL;
         }
        close POT;
        system ("rm -f POTCAR.*");
}
#######################################################"

&main ();
