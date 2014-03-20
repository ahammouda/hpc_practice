#!/bin/bash
touch global_out2d.dat
for (( c=0; c<10; c++ ))
do
   echo "Welcome $c times"
   cat $c.out2d.dat >> global_out2d.dat
done
