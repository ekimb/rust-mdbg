cat $1 | grep -v target | grep -o -E "\"utg[^']*\""  | tr -d \" |sort|uniq > unitigs.$1.txt
grep -f unitigs.$1.txt unitigs_genome_assoc.txt > unitigs.$1.data.txt
python get_genomes_ids.py unitigs.$1.data.txt > unitigs.$1.genomes.txt
echo "type in a internet-friendly machine:"
echo "pysradb metadata --saveto unitigs.$1.ena.txt  \$(< unitigs.$1.genomes.txt)"
read -p "press any key"
cut -f 6 unitigs.$1.ena.txt
echo "distinct taxids: $(cut -f 5 unitigs.$1.ena.txt |sort|uniq|wc -l)"

