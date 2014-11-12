#!/bin/bash
#===========================================================
# 1) extract cv1.txt including multicomponent data
# 2) plot extracted data as xy graph
# 3) merge extracted data for 3D plotting
# 4) plot phasediagram
#
#
# usage; cvmplt.sh
#
# 18/Mar/2011
# 18/Mar/2012 reconstruct for mainly multicomponent  
# 29/Mar/2012 modified for volume
#===========================================================

GNUPLOT="gnuplot"
NA=6.02e23
kB=1.380662e-23

#-----------------------------------------------------------
# check cv1.txt
# add description for each row
#-----------------------------------------------------------

count=0

while read cur_line;
do   
   set -- $cur_line
   prefix=$1

   case $prefix
   in
   f)
   count=`expr $count + 1`
   if [ $count = 2 ]; then
      break   
   fi
   continue
    ;;  

   c)
   row=$#
   row=`expr $row - 1`
   if [ $row -ge 2 ]; then
      echo ""
      echo "-----------------------------------------------------------"
      echo " This is the multi-component system !!"
      echo " Number of independent composition is $row"
      echo "-----------------------------------------------------------"
      echo ""
   fi
   tmpfilelist1="$tmpfilelist1 cv1/p\${num_phase1}_t\${temp1}c"
   tmpfilelist1="$tmpfilelist1 cv1/p\${num_phase1}_t\${temp1}f"
   tmpfilelist="$tmpfilelist cv1/p\${num_phase}_t\${temp}c"
   tmpfilelist="$tmpfilelist cv1/p\${num_phase}_t\${temp}f"
    ;;

   m)
   tmpfilelist1="$tmpfilelist1 cv1/p\${num_phase1}_t\${temp1}m"
   tmpfilelist="$tmpfilelist cv1/p\${num_phase}_t\${temp}m"
   continue
    ;;
     
   v)
   echo ""
   echo "volume was found"
   echo ""
   tmpfilelist1="$tmpfilelist1 cv1/p\${num_phase1}_t\${temp1}v"
   tmpfilelist="$tmpfilelist cv1/p\${num_phase}_t\${temp}v"
    ;;

   e)
   tmpfilelist1="$tmpfilelist1 cv1/p\${num_phase1}_t\${temp1}e"
   tmpfilelist="$tmpfilelist cv1/p\${num_phase}_t\${temp}e"
   continue
    ;;
   esac
done < cv1.txt

#echo "$tmpfilelist"

#-----------------------------------------------------------
# start main menu
#-----------------------------------------------------------
menu="ext-cv1 setdata xy-plot merge F-3Dplot Tx-plot exit"
select command in $menu
do

  if [ "$command" = "ext-cv1" ]; then

  if [ -d cv1 ]; then
#     echo " " > /dev/null
     rm -rf cv1
     mkdir cv1
  else
     mkdir cv1
  fi
  
  while read cur_line;
  do 
   
    set -- $cur_line

    prefix=$1

    case $prefix
    in
 
     f)
     # store num_phase & temp last cycle
     num_phase1=$num_phase
     temp1=$temp

     # store num_phase & temp current cycle
     num_phase=$2
     temp=$3
     OP1=$4
     OP2=$5
     F=$6
     S=$7
     cond=$8
     num_corr=$9
     num_prob=${10}
     unknown=${11}
     sum_deriv=${12}
     
     # writing output file about 1 cycle before
     if [ -e cv1/p"${num_phase1}"_t"${temp1}"f ]; then
     filelist1=`eval echo $"$tmpfilelist1"`
     paste -d " " $filelist1 >cv1/p"${num_phase1}"_t"${temp1}"
     fi

     # writing output temp file about this cycle
     if [ -e cv1/p"${num_phase}"_t"${temp}"f ]; then
        echo $OP1 $OP2 $F $S $cond $num_corr $num_prob $unknown $sum_deriv >> cv1/p"${num_phase}"_t"${temp}"f
     else
        echo "OP1 OP2 F S cond num_corr num_prob unknown sum_deriv " >> cv1/p"${num_phase}"_t"${temp}"f
        echo $OP1 $OP2 $F $S $cond $num_corr $num_prob $unknown $sum_deriv >> cv1/p"${num_phase}"_t"${temp}"f
     fi
      ;;  

     c)
     head_comp=""
     line_comp=""
     i=2
     while [ ${i} -le $# ]
     do
        num_comp=`expr ${i} - 1`
        comp[$num_comp]=`eval echo '$'${i}`
#        echo ${comp[1]} ${comp[2]}
        head_comp="$head_comp comp${num_comp}"
        line_comp="$line_comp ${comp[${num_comp}]}"
        i=`expr ${i} + 1`
     done

     if [ -e cv1/p"${num_phase}"_t"${temp}"c ]; then
        echo $line_comp >> cv1/p"${num_phase}"_t"${temp}"c
     else
        echo "#$head_comp" >> cv1/p"${num_phase}"_t"${temp}"c
        echo $line_comp >> cv1/p"${num_phase}"_t"${temp}"c
     fi
      ;;

     m)
     chem_pot=$2
     if [ -e cv1/p"${num_phase}"_t"${temp}"m ]; then
        echo $chem_pot >> cv1/p"${num_phase}"_t"${temp}"m
     else
        echo "chem_pot " >> cv1/p"${num_phase}"_t"${temp}"m
        echo $chem_pot >> cv1/p"${num_phase}"_t"${temp}"m
     fi
      ;;
     
     v)
     vol=$2
     if [ -e cv1/p"${num_phase}"_t"${temp}"v ]; then
        echo $vol >> cv1/p"${num_phase}"_t"${temp}"v
     else
        echo "volume " >> cv1/p"${num_phase}"_t"${temp}"v
        echo $vol >> cv1/p"${num_phase}"_t"${temp}"v
     fi
      ;;
     e)
     eigen=$2
     if [ -e cv1/p"${num_phase}"_t"${temp}"e ]; then
        echo $eigen >> cv1/p"${num_phase}"_t"${temp}"e
     else
        echo "eigen" >> cv1/p"${num_phase}"_t"${temp}"e
        echo $eigen >> cv1/p"${num_phase}"_t"${temp}"e
     fi
      ;;
   esac

  done < cv1.txt

  filelist=`eval echo $"$tmpfilelist"`
  paste -d " " $filelist >cv1/p"${num_phase}"_t"${temp}"
  rm -rf cv1/*.c cv1/*.f cv1/*.m cv1/*.e cv1/*.v
  
  elif [ "$command" = "setdata" ]; then



  cd cv1

  #-----------------------------------------------------------
  # extract menu from p#_t# file
  #-----------------------------------------------------------

  if [ -f cv1plot.gnu ]; then
   rm cv1plot.gnu
  fi

  files=`ls p*`
  set -- $files
  menu=`head -n 1 $1`

  echo ""
  echo "# available items"
  echo "$menu"
  
  echo "# input row number for x axis"
  read index_x
  
  echo "# input row number for y axis"
  read index_y
  
  echo "# need to change the scale from K to kJ/mol ? (y/n)"
  read index_change
  if [ $index_change = y ]; then
     index_y=\(\$${index_y}*$kB*$NA/1000\)
  fi

  echo ""
  echo "# input phase"
  read phase

  echo "# input temperatures (all=>a)"
  read temp_plot
  set -- $temp_plot

  if [ "$1" = "a" ];then

    list=(`ls | grep "p${phase}_"`)

    for cur_file in "${list[@]}"
    do
      file_list=$file_list", \"$cur_file\" u $index_x:$index_y"

    done
    
  else
    while [ "$1" != "" ]
    do

      list=(`ls | grep "p${phase}_" | grep "t${1}\."`)

      for cur_file in "${list[@]}"
      do
        file_list=$file_list", \"$cur_file\" u $index_x:$index_y"

      done
      shift
    done
  fi

  echo $file_list

  cd ..

  elif [ "$command" = "xy-plot" ]; then


  cd cv1

########################################################################
cat - > cv1plot.gnu <<END
set terminal x11

set xlabel ""
set ylabel ""
set title ""
set pointsize 2
set xzeroaxis
set xtics autofreq
plot 0$file_list
pause -1
END
########################################################################

  $GNUPLOT cv1plot.gnu

  cd ..

  elif [ "$command" = "merge" ]; then

  if [ -f merge.dat ]; then
   rm merge.dat
   echo "data.dat is deleted"
  fi

  availdir=`ls -F | grep /`
  echo ""
  echo "# available items"
  echo $availdir

  echo "# input dirs to merge (note; space delimiter & including "/")"
  read mergedirs
  set -- $mergedirs
  
  echo "# input phase"
  read phase

  echo "# input temperature"
  read temp

  echo "" >space
#  file_list=""
  while [ "$1" != "" ]
  do
#     file_list="$file_list $1p${phase}_t${temp}."
     cat $1p${phase}_t${temp}. space >>merge.dat
     shift
  done

#  cat $file_list > data.dat
   rm space

  elif [ "$command" = "F-3Dplot" ]; then

########################################################################
cat - > ternary.gnu <<END
#set terminal postscript color
set terminal X11
set view 66,5
set bmargin 3
set lmargin 3
set rmargin 3
set tmargin 3
set size ratio 0.866
set yrange [0:0.866]
set xrange [0:1]
set noborder
set noxtics
set noytics
set label at graph -0.1,-0.1,0.5 "kJ/mol" rotate by 90
set label '' at 0, -0.03,10 center
set label '' at 1, -0.03,10 center
set label '' at 0.5, 0.886,10 center
set output 'ternary.eps'

set style line 1 lt 1 lw 3 pt -1 ps 1
set style line 2 lt 5 lw 1 pt -1 ps 1

# x
set arrow 1 from 0.0,0,-35 to 1.00, 0.000,-35 nohead
set arrow 2 from 0.1,0,-35 to 0.55, 0.779,-35 nohead
set arrow 3 from 0.2,0,-35 to 0.60, 0.693,-35 nohead
set arrow 4 from 0.3,0,-35 to 0.65, 0.606,-35 nohead
set arrow 5 from 0.4,0,-35 to 0.70, 0.520,-35 nohead
set arrow 6 from 0.5,0,-35 to 0.75, 0.433,-35 nohead
set arrow 7 from 0.6,0,-35 to 0.80, 0.346,-35 nohead
set arrow 8 from 0.7,0,-35 to 0.85, 0.260,-35 nohead
set arrow 9 from 0.8,0,-35 to 0.90, 0.173,-35 nohead
set arrow 10 from 0.9,0,-35 to 0.95, 0.0866,-35 nohead

# z
set arrow 11 from 1.00, 0.0000,-35 to 0.50, 0.8660,-35 nohead
set arrow 12 from 0.95, 0.0866,-35 to 0.05, 0.0866,-35 nohead
set arrow 13 from 0.90, 0.1730,-35 to 0.10, 0.1730,-35 nohead
set arrow 14 from 0.85, 0.2600,-35 to 0.15, 0.2600,-35 nohead
set arrow 15 from 0.80, 0.3460,-35 to 0.20, 0.3460,-35 nohead
set arrow 16 from 0.75, 0.4330,-35 to 0.25, 0.4330,-35 nohead
set arrow 17 from 0.70, 0.5200,-35 to 0.30, 0.5200,-35 nohead
set arrow 18 from 0.65, 0.6060,-35 to 0.35, 0.6060,-35 nohead
set arrow 19 from 0.60, 0.6930,-35 to 0.40, 0.6930,-35 nohead
set arrow 20 from 0.55, 0.7790,-35 to 0.45, 0.7790,-35 nohead

# y
set arrow 21 from 0.50, 0.866 ,-35 to 0.0,0,-35 nohead
set arrow 22 from 0.05, 0.0866,-35 to 0.1,0,-35 nohead
set arrow 23 from 0.10, 0.173 ,-35 to 0.2,0,-35 nohead
set arrow 24 from 0.15, 0.260 ,-35 to 0.3,0,-35 nohead
set arrow 25 from 0.20, 0.346 ,-35 to 0.4,0,-35 nohead
set arrow 26 from 0.25, 0.433 ,-35 to 0.5,0,-35 nohead
set arrow 27 from 0.30, 0.520 ,-35 to 0.6,0,-35 nohead
set arrow 28 from 0.35, 0.606 ,-35 to 0.7,0,-35 nohead
set arrow 29 from 0.40, 0.693 ,-35 to 0.8,0,-35 nohead
set arrow 30 from 0.45, 0.779 ,-35 to 0.9,0,-35 nohead

# for frame
set arrow 31 from 0,0,-35 to 0,0,5 nohead
set arrow 32 from 0.5,0.866,-35 to 0.5,0.8660,5 nohead
set arrow 33 from 1.0,0,-35 to 1.0,0,5 nohead
set arrow 34 from 0,0,5 to 1, 0.0,5 nohead
set arrow 35 from 1, 0,5 to 0.50, 0.866,5 nohead
set arrow 36 from 0.50, 0.866,5 to 0,0,5 nohead


NA=6.02e23
kB=1.380662e-23

set parametric
set pm3d
END
echo 'splot "merge.dat" u ($1+2*(1-$1-$2))/(2*($1+$2+(1-$1-$2))):(sqrt(3)*$1/(2*($1+$2+(1-$1-$2)))):($4*kB*NA/1000) notitle with lines' >> ternary.gnu
echo 'pause -1' >>ternary.gnu
########################################################################

  $GNUPLOT ternary.gnu

  elif [ "$command" = "Tx-plot" ]; then
########################################################################
cat -> phasediagram.gnu << END
#!/opt/local/bin/gnuplot -persist
set terminal x11
#set terminal postscript color
#set output "phasediagram.eps"

set pointsize 1
set title "Phase Diagram" 
set xlabel "x" 
set xrange [ 0.00000 : 1.00000 ] 
set ylabel "T (K)" 
#set yrange [ 0.00000 : 1000.00 ] 
set locale "ja_JP.UTF-8"
plot "1phaseregion.txt" u 2:1, "2phaseregion.txt" u 2:1, "2phaseregion.txt" u 3:1, "spinodalregion.txt" u 2:1
pause -1
#    EOF
END
########################################################################

   $GNUPLOT phasediagram.gnu

elif [ "$command" = "exit" ]; then
    clear
    exit

else
    clear
    echo $REPLY is not in the command list.
    echo
fi
done





