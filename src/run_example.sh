echo -e ""
echo -e "\t+----------"
echo -e "\t| Jordi Lladós Segura"
echo -e "\t| Memory-Efficient Consistency Library for Multiple Sequence Alignment Tools"
echo -e "\t+----------"
echo -e ""
echo -e "this script launch t_coffee with additional parameters."

#INPUT
SEQ=../example/GEL_100
TREE=../TREES/GEL_100.dnd

#PARAMETERS
EX=very_fast_triplet
DP=myers_miller_pair_wise

./t_coffee $SEQ -usetree=$TREE -extend_mode=$EX -dp_mode=$DP
