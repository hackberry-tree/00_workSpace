#!/usr/bin/perl
# mixing energies from SQS calculation
# ver_1.91 - may 2013

use warnings;

# where are pure files ?
$michi='/home/enoki/zengen2.0/pure/';                              
#ex: $michi='/home/jcc/bin/pure/';  

# Variables
        $nA = 16;         # Nb atoms by cell
        $nE = 2;          # Nb of differents elements
        $ECUT = 400;      # Cutoff-energy

##########################################################
sub main
{
        my @ld = glob ('*');
        @compo = ();      # table of compositions studied
        @energ = ();      # table of total energy         
        @vol = ();        # table of volume               
        @El = ();         # table of elements name
        @SER = ();        # table of SER energie
        @SRE = ();        # table of SRE energie
#        @HFE = ();        # table of Heat-formation Energie
#        @ME = ();         # table of Mixing Energie

        system("clear");
        my @states = ("vol", "full");   # vol full  MAG ?
        foreach $STATE(@states) 
        {
	 print "#######################################\n";
         print "#  relaxation step: $STATE\n";
         foreach $phase (@ld)
         {
         if (-d $phase)
           {
	   print "#######################################\n";
           print "$phase  - Ecut=$ECUT \n";
	   print "#######################################\n";
           chdir "$phase";
                $nC = 0;        # Number of differents compositions in $dir
                &looksum ();    # Find the energy and volume of each composition
                &findelem ();   # Find the name of each element         
	  print "#######################################\n";
          my $loulou = 0;       # temporary variable  
          while ($loulou < $nE) # For each element, find the SER and SRE of $phase  
          {   
           $XB=100*$loulou;
           $AT=$El[$loulou];                                  # pour l'element AT
           $SER[$loulou]=&lookSER ;                           # energie SER est : 
           $SRE[$loulou]=&lookSRE ;                           # energie SRE est : 
           $loulou += 1;
	  print "#######################################\n";
          }
          &CalcHeat ();
	  print "#######################################\n";
          chdir "..";
           }
         }
        }
}


##########################################################
sub looksum
{
                open (SUM,"sum.$STATE.out");
                while (<SUM>)
                {
                        if ($_ =~ m/:/)
                        {
                                my @ligne = split ();
                                $compo[$nC] = $ligne[0];
                                $energ[$nC] = $ligne[2];
                                $vol[$nC]   = $ligne[7];
                                $nC++;
                        }
                }
                close SUM;
}
##########################################################
sub findelem
{
                $titi = 0;
                open (POT,"POTCAR");
                while (<POT>)
                {
                        if ($_ =~ /^\s*(PAW_PBE)/)  #^commence par \s* même * espace
                        {
                                my @ligne = split ();
                                $El[$titi]  = $ligne[1];
                                $titi++;
                        }
                }
                close POT;
}
##########################################################
# looking for the Stable Element Reference energy ########
sub lookSER   
{
               open (SER,"$michi/SER-$ECUT.out");
               while (<SER>)
                {
                        if ($_ =~ m/$AT/)
                        {
                                my @ligne = split ();
                                $SS = $ligne[1];
                                $SE = $ligne[3];
                                $SM = $ligne[6];
                                $SV = $ligne[8];
                        }
                }
                close SER;
                print " $AT - $SS - SER,     E=$SE eV, vol= $SV A3/at, M=$SM\n";
                $SE;
}
##########################################################
# looking for the Structure Refered Energy of phase ######
sub lookSRE
{
               open (SER,"$michi/$phase-$ECUT.out");
               while (<SER>)
                {
                        if ($_ =~ m/$AT/)
                        {
                                my @ligne = split ();
                                $SE = $ligne[2];
                                $SV = $ligne[7];
                        }
                }
                close SER;
                print " $AT - $phase,           E=$SE eV, vol= $SV A3/at \n";
                $compo[$nC] = $XB ; $energ[$nC] = 16*$SE; $vol[$nC] = 16*$SV ; $nC++; # add uniaries in Table
                $SE;
}
##########################################################
# Calculation of heat of formation and heat of mixing
sub CalcHeat
{
           $toto=0;
           system("rm -f ../$El[0]-$El[1].$phase.$STATE.dat");
           open (DAT,">> ../$El[0]-$El[1].$phase.$STATE.dat");
           print "# x_$El[1]   F(eV)    |         Heat-formation      |      Heat-mixing            | volume (A^3/at)\n" ;
           print DAT "# x_$El[1]   F(eV)    |         Heat-formation      |      Heat-mixing            | volume (A^3/at)\n" ;
           while ($toto < $nC)
           {
            #print $energy[$toto];
            $E = $energ[$toto] / $nA;
            $H = $E - (($compo[$toto])/100) * $SER[1] ;
            $H = $H - (1-($compo[$toto])/100) * $SER[0] ;
            $HKJ = $H*96.486;
            $M = $E - (($compo[$toto])/100) * $SRE[1] ;
            $M = $M - (1-($compo[$toto])/100) * $SRE[0] ;
            $MKJ = $M*96.486;
#            print "$toto $compo[$toto] $energ[$toto] | $H eV $HKJ kJ/mole | $M eV $MKJ kJ/mole\n" ;
            printf "%4d %12.5f | %8.4f eV %8.4f kJ/mol | %8.4f eV %8.4f kJ/mol | %8.4f\n", $compo[$toto], $energ[$toto], $H, $HKJ, $M, $MKJ, $vol[$toto]/$nA  ;
            printf DAT "%4d %12.5f | %8.4f eV %8.4f kJ/mol | %8.4f eV %8.4f kJ/mol | %8.4f\n", $compo[$toto], $energ[$toto], $H, $HKJ, $M, $MKJ, $vol[$toto]/$nA  ;
            $toto += 1;
           }
            close DAT;
}


&main ();
