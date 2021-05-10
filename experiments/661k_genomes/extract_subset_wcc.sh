n=1000 # too much for gephi
n=10
echo "<?xml version='1.0' encoding='utf-8'?>
<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">
  <graph edgedefault=\"directed\">"   > wcc_subset.graphml
head -n $n reggraph-k10-p0.001-l12.u_wcc.sizes  |awk '{print $2}' |xargs cat |grep -v key |grep -v graph |grep -v xml >> wcc_subset.graphml
echo  "  </graph>
</graphml>" >> wcc_subset.graphml
