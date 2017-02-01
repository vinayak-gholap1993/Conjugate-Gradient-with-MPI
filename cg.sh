#!/bin/bash

make clear
clear
make

echo

read -p "Number of Processors : " numProc
read -p "Number of grids in X-direction : " nx
read -p "Number of grids in Y-direction : " ny
read -p "Number of iterations : " iter

mpirun -np $numProc ./cg $nx $ny $iter 

echo 
