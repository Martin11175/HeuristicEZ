#!/bin/bash
# 
# Execute a complete array of proximity tests in matlab
# Author: mh881@york.ac.uk

# Algorithm part inclusion
rgea_types=( "'ground'" "'rgea'" "'simple'" "'none'" )
ap_types=( "'all'" "'overlap'" "'max'" "'none'" )
loc_types=( "'all'" "'overlap'" "'max'" "'none'" )
ips_types=( "'all'" ) #( "'all'" "'max'" )

# Model generation parameters
thresholds="[-100 -95 -90 -85 -80 -75 -70 -65 -60 -55 -50]"
buildings="[0 1 2]"
floors="[0 1 2 3]"
bounds=( 20 ) #( 0 5 10 15 20 25 30 40 50 75 100 )
chances=( 100 ) #( 100 75 50 30 25 20 15 10 5 )

# Do every combination of parameters and parts
for rgea_type in ${rgea_types[@]}
do
    for ap_type in ${ap_types[@]}
    do
        for loc_type in ${loc_types[@]}
        do
            for ips_type in ${ips_types[@]}
            do
                for bound in "${bounds[@]}"
                do
                    for chance in "${chances[@]}"
                    do
                        echo "Running: " $rgea_type $ap_type $loc_type $ips_type $bound $chance
                        # Run tests
                        nohup matlab -nosplash -nojvm -nodisplay -r "\
                        rgea_type = $rgea_type;\
                        ap_type = $ap_type;\
                        loc_type = $loc_type;\
                        ips_type = $ips_type;\
                        bounds = $bound;\
                        chance = $chance;\
                        buildings = $buildings;\
                        floors = $floors;\
                        thresholds = $thresholds;\
                        main;\
                        exit" >/dev/null 2>&1 &
                    done
                done
            done
        done
    done
done
