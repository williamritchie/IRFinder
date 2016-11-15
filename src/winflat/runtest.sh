#!/bin/sh
echo -e " " \\t 2E=0.01 \\t 2E=0.005
echo -e x \\t Ymin--Ymax \\t Ymin--Ymax
for i in 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20
do
echo -e  $i \\t `./winflat -xvalue $i -sig 0.01 | awk '{print $2}'` \\t\\t  `./winflat -xvalue $i -sig 0.05 | awk '{print $2}'` 
done
exit 0 
