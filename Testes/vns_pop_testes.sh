#!/bin/bash

# codigos de funcoes do VNS
FUNCTIONS_CODE=("0" "1" "2" "3" "4" "5" "7" "12" "13" "77")
FUNCTIONS=( "RASTRIGIN" "SCHAFFER_F7" "GRIEWANK" "ACKLEY" "ROSENBROCK" "SPHERE" "SCHAFFER_F6" "LEVY" "ZAKHAROV"  "SCHWEFEL222" )

DIMENSIONS=("250")

METHODS_CODE=( "6" )
METHODS=( "PRVNS" )

P_CROSS=("0.9" )

	echo "Creating dirs ..."
	echo "Generating input file ..."

	k=0;
	for p in "${P_CROSS[@]}"
	do
		j=0;
		for m in "${METHODS[@]}"
		do
			i=0;
			for f in "${FUNCTIONS[@]}"
			do
				for d in "${DIMENSIONS[@]}"
				do
					mkdir -p $m/$f/$d
					cp input.in $m/$f/$d/
					cd $m/$f/$d/
					sed -i.bak s/FUNCTIONS/${FUNCTIONS_CODE[$i]}/g input.in
					rm input.in.bak
					sed -i.bak s/METHODS/${METHODS_CODE[$j]}/g input.in
					rm input.in.bak
					sed -i.bak s/P_CROSS/${P_CROSS[$k]}/g input.in
					rm input.in.bak
					cd ../../../
				done
				let i++;
			done
			let j++;
		done
		let k++;
	done
	
	echo "Generation complete!"
	echo "Runing all ..."
	for m in "${METHODS[@]}"
	do
		for f in "${FUNCTIONS[@]}"
		do
			# mkdir Convergencia
			#wait finish before run more testes
			while pgrep "alg" > /dev/null; do sleep 10; done
				gnome-terminal -x sh -c "{ time ./alg $m/$f/250/input.in > $m/$f/250/output.out ; } 2>> $m/$f/250/output.out ;"
				sleep 1
			while pgrep "alg" > /dev/null; do sleep 10; done

			# media=`ls Convergencia/*Media.data`
			# melhor=`ls Convergencia/*Melhor.data`
			# saida_img=$m-$f-250.png
			# ./gnuplot_grafico.sh $saida_img $media $melhor
			# mv $saida_img Convergencia
			# mv Convergencia $m/$f/250/
		done
	done
