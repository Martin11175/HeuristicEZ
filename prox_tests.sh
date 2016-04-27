#!/bin/bash
# 
# Execute a complete array of proximity tests in matlab
# Author: mh881@york.ac.uk

# EZ params
ez_prox_threshold=( 1 2 3 4 5 6 7 8 ) # Proximity threshold in dB
min_AP_overlap=( 2 3 4 5 6 7 ) # Minimum similar visible APs

# Heuristic params
min_AP_strength=( -100 -90 -80 ) #( -100 -95 -90 -85 -80 -75 -70 -65 -60 )
num_strongest=( 2 3 4 5 7 10 ) #( 2 3 4 5 6 7 8 9 10 )
min_strong_overlap=( 2 3 4 5 6 ) #( 1 2 3 4 5 6 7 8 9 10 )

# Parameter counters
ez_prox_threshold_n=0
min_AP_overlap_n=0
min_AP_strength_n=0
num_strongest_n=0
min_strong_overlap_n=0

# Iterate until all parameters for each metric are tested
more_ex_tests=true
more_h_test=true

while $more_ex_tests || $more_h_tests
do
    echo "Starting tests: (" $ez_prox_threshold_n $min_AP_overlap_n ") ("\
        $min_AP_strength_n $num_strongest_n $min_strong_overlap_n ")"
    # Run tests
    #nice -n -5 \
    nohup \
    matlab -nosplash -nojvm -nodisplay -r "\
    ez_prox_threshold = ${ez_prox_threshold[$ez_prox_threshold_n]};\
    min_AP_overlap = ${min_AP_overlap[$min_AP_overlap_n]};\
    min_AP_strength = ${min_AP_strength[$min_AP_strength_n]};\
    num_strongest = ${num_strongest[$num_strongest_n]};\
    min_strong_overlap = ${min_strong_overlap[$min_strong_overlap_n]};\
    prox_test;\
    exit" >/dev/null 2>&1 &

    # Increment parameters
    if $more_ex_tests
    then
        ez_prox_threshold_n=$(($ez_prox_threshold_n + 1))
        if [ $ez_prox_threshold_n -eq ${#ez_prox_threshold[@]} ]
        then
            ez_prox_threshold_n=0
            min_AP_overlap_n=$(($min_AP_overlap_n + 1))
        fi

        if [ $min_AP_overlap_n -eq ${#min_AP_overlap[@]} ]
        then
            more_ex_tests=false
            min_AP_overlap_n=0
        fi
    fi

    if $more_h_tests
    then
        min_AP_strength_n=$(($min_AP_strength_n + 1))
        if [ $min_AP_strength_n -eq ${#min_AP_strength[@]} ]
        then
            min_AP_strength_n=0
            num_strongest_n=$(($num_strongest_n + 1))
        fi

        if [ $num_strongest_n -eq ${#num_strongest[@]} ]
        then
            num_strongest_n=0
            min_strong_overlap_n=$(($min_strong_overlap_n + 1))
        fi

        if [ $min_strong_overlap_n -eq ${#min_strong_overlap[@]} ]
        then
            more_h_tests=false
            min_strong_overlap_n=0
        fi
    fi
done
