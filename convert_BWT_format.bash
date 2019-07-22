input=$1
output=$2

cat $input/bwt_test.0.dat | tr "\000\001\002\003\004\005" "\$ACGNT" | msbwt convert $output

