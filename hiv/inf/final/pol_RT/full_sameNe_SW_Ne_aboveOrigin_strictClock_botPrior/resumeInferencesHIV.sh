module load java
for i in {0..2}; do
bsub -J "hiv_bp" -W 24:00 -R "rusage[mem=10000]"  java -Xmx42g -jar ../../TnT.jar -seed random -resume env_pol_saConstrained_s1_sameNe_run${i}.xml
done
