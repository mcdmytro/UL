> results.txt
> results_short.txt

declare -A comb
comb[0,0]=2317
comb[0,1]=4274
comb[1,0]=2317
comb[1,1]=4685
comb[2,0]=2460
comb[2,1]=4630
comb[3,0]=2460
comb[3,1]=4500
comb[4,0]=2460
comb[4,1]=4700

declare -A varns
varns[0,0]=0
varns[0,1]=0
varns[1,0]=1
varns[1,1]=0
varns[2,0]=-1
varns[2,1]=0
varns[3,0]=0
varns[3,1]=1
varns[4,0]=0
varns[4,1]=-1

for a in 0 1 2 3 4
do
    for b in 0 1 2 3 4 
    do
        > config.h

        echo "#if !defined(MYLIB_CONSTANTS_H)" >> config.h
        echo "#define MYLIB_CONSTANTS_H 1" >> config.h
        echo " " >> config.h
        echo "int decay = ${comb[$a,0]};" >> config.h
        echo "int state = ${comb[$a,1]};" >> config.h
        echo " " >> config.h
        echo "// For systematics evaluation" >> config.h
        echo "double mass_err_fact  = ${varns[$b,0]};" >> config.h
        echo "double width_err_fact = ${varns[$b,1]};" >> config.h
        echo "" >> config.h
        echo "#endif" >> config.h
        # echo "${varns[$b,0]} ${varns[$b,1]}"

        ./run_UL_calc.sh
    done
done

root -l errors_calc.C