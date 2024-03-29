#!/bin/bash

#0.- Main path
mainpath="../../2_equilibration"
          
#1.- Directory where the data is
datapath=${mainpath}"/data/"

#2.- Directory where the output will be stored
mkdir -p ${mainpath}"/cv_vs_time"

outputpath=${mainpath}"/cv_vs_time/"

#3.- Are we reading multiple configurations? (1-yes)
multi=0

#4.- which configuration do you want to read?
conf=4

#5.- Starting timestep
start=0

#6.- The dump frequency
dumpfreq=100000

#7.- Last timestep of the simulation
last=400000000

#8.- number Nucluotides in the system
N=147

#9.- oxDNA version
version=2


#####################
#compile the program#
#####################
c++ cv_evolution.cpp -o cv_evolution.out -lm



#################
#RUN THE PROGRAMM#
##################
./cv_evolution.out ${datapath} ${outputpath} ${multi} ${conf} ${start} ${dumpfreq} ${last} ${N} ${version}
