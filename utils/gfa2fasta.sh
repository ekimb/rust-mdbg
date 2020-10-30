awk '/^S/{print ">"$2"\n"$3}' $1.gfa | fold > $1.fa
