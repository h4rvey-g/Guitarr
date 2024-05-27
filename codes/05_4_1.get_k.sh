source miniQuant/base/bin/activate
python miniQuant/isoform_quantification/main.py cal_K_value \
  -gtf data/reference/GRCh38_chr.gtf \
  -t 60 \
  -o data/reference/GRCh38_chr_K_value
# python miniQuant/isoform_quantification/main.py cal_K_value \
#   -gtf "data/reference/3rd_gen/3rd_gen.gtf" \
#   -t 60 \
#   -o data/reference/3rd_gen_K_value
# prefixes=("GS689" "PC3E")
# suffixes=("_1" "_2" "_3")
# combinations=()
# for prefix in "${prefixes[@]}"; do
#   for suffix in "${suffixes[@]}"; do
#     combinations+=("${prefix}${suffix}")
#   done
# done

# root=data/04.Tx_quantification/01.Stringtie/
# for combination in "${combinations[@]}"; do
#   python miniQuant/isoform_quantification/main.py cal_K_value \
#     -gtf ${root}/${combination}/${combination}.gtf \
#     -t 40 \
#     -o ${root}/${combination}/${combination}_K_value.txt
# done


# root=data/04.Tx_quantification/02.Stringtie_with_impute/
# # generate combinations of GS689, PC3E and _1, _2, _3
# for combination in "${combinations[@]}"; do
#   python miniQuant/isoform_quantification/main.py cal_K_value \
#     -gtf ${root}/${combination}Aligned.sortedByCoord.out.imputed.bam/${combination}Aligned.sortedByCoord.out.imputed.bam.gtf \
#     -t 40 \
#     -o ${root}/${combination}Aligned.sortedByCoord.out.imputed.bam/${combination}_K_value.txt
# done

# root=data/04.Tx_quantification/02.Stringtie_with_impute/bam_impute

# for combination in "${combinations[@]}"; do
#   python miniQuant/isoform_quantification/main.py cal_K_value \
#     -gtf ${root}/${combination}_Aligned.sortedByCoord.out/
#     -t 40 \
#     -o ${root}/${combination}/${combination}_K_value.txt
# done