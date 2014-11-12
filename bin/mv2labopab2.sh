#!/bin/bash
path=/labopub2`pwd | sed 's/\/home//'`
fname=`echo $1 | sed 's/\///'`.tar.gz
dst=${path}/${fname}
mkdir -p $path
tar cvfz $dst $1

rm -r $1/*
touch $1/mv2labopub2
