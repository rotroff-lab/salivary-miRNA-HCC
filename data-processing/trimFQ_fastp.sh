#!/bin/sh
#SBATCH --job-name="fastp"
#SBATCH --error=slurm/fastp-sbatch-%j.error
#SBATCH --output=slurm/fastp-sbatch-%j.out


module load trim_galore
module load fastp


samples=("00_913A_14_S1_R1" "01_656B_22_S2_R1" "09_497B_31_S3_R1" "11_902A_10_S4_R1" \ 
         "23_918A_17_S5_R1" "28_518B_25_S6_R1" "28_921A_19_S7_R1" "30_585B_29_S8_R1" \ 
         "36_671B_33_S9_R1" "38_917A_16_S10_R1" "40_934B_28_S11_R1" "44_613B_36_S12_R1" "50_490B_21_S13_R1" \ 
         "50_887A_02_S14_R1" "55_517B_34_S15_R1" "55_644B_26_S16_R1" "58_910A_12_S17_R1" \ 
         "59_599B_35_S18_R1" "59_666B_38_S19_R1" "60_760B_30_S20_R1" "64_899A_07_S21_R1" \ 
         "68_912A_13_S22_R1" "71_664B_23_S23_R1" "72_635B_27_S24_R1" "73_747B_40_S25_R1" \ 
         "75_652B_37_S26_R1" "77_708B_39_S27_R1" "77_901A_09_S28_R1" "77_922A_20_S29_R1" \ 
         "79_904A_11_S30_R1" "80_923A_05_repeat_S31_R1" "82_920A_18_S32_R1" "84_928A_04_S33_R1" \ 
         "87_888A_03_repeat_S34_R1" "89_555B_24_S35_R1" "92_740A_01_repeat_S36_R1" \ 
         "93_898A_06_S37_R1" "93_900A_08_S38_R1" "98_509B_32_S39_R1" "99_914A_15_S40_R1")


for s in ${samples[@]}
do
	echo "$s"
        fastp --in1 042220_NextSeq/${s}.fastq.gz \
        --out1 trimmed/${s}_TRIMMED.fastq.gz \
        -L -q 20 -h "${s}_fastpReport.html" \
        --overrepresentation_analysis 
done


