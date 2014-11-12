#!/bin/bash
# Execute RIETAN
# The environment variable RIETAN should be defined before execution
f=${1##*/}
f=${f%.ins}
SAMPLE_DIR=${1%/*}
RIETAN=${2%/*}
export RIETAN
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
export LD_LIBRARY_PATH
cd "$SAMPLE_DIR"
$2 $f.ins $f.int $f.bkg $f.itx $f.hkl $f.xyz $f.fos $f.ffe $f.fba $f.ffi $f.ffo $f.vesta $f.plt $f.gpd > $f.lst

if [ -e $f.itx ]
 then
	if [ -z $3 ] || [ "$3" = "\"" ]
	 then
		xdg-open $f.itx
	else
		dir=${3%/*}
		ext=${3##*.}
		if [ "$ext" = "jar" ]
		 then
			exe=${3##*/}
			cd $dir
			java -jar $exe "$SAMPLE_DIR/$f.itx" &
		else
			$3 $f.itx
		fi
	fi
else
	xdg-open "$SAMPLE_DIR/$f.lst"
fi

#sleep 5
#rm $SAMPLE_DIR/$f.*
