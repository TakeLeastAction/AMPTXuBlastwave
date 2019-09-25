#!/bin/bash
rm uni_ball_small.dat
rm NJLMUL.dat
touch uni_ball_small.dat 
rm ampt
rm zpcBW.dat
ifort -o ball ball_ini.f90
#f77 -o ampt -O *.f -mcmodel=medium

f77 -o ampt -O main.f zpc.f NJL.f linana.f hipyset1.35.f hijing1.383_ampt.f art1f.f amptsub.f -mcmodel=medium

