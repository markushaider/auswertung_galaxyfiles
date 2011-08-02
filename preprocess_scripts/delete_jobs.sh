#!/bin/bash

# Start Number of job
start=28598
end=30000

while test $start -le $end; do
	qdel $start
	wait
	start=$(( $start + 1 ))
done

