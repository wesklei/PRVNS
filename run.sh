#!/bin/bash
# gcc -c mersenne.c -lm
# gcc mersenne.o -o alg algorithm.c -lm
make clean
make
mv alg Testes
cd Testes
chmod +x vns_pop_testes.sh
./vns_pop_testes.sh
