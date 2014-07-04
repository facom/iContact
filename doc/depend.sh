for theta in $(seq 1 2 90)
do
    echo -n "$theta "
    python elipsoid.py $theta
done
