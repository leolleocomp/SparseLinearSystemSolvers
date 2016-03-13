#!/bin/bash

for file in input/*
do
	name=${file##*/}
	name=${name%%.*}

	./novos_testes "resultados/${name}_iterações.out"\
		       "resultados/${name}_tempo.out" < "$file" 2> "resultados/log/${name}" 
done
