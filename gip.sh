#!/bin/bash

# Default values
num_cpus=0
num_gpus=0
ref_panel=""

# Function to display usage
usage() {
    echo "Usage: $0 -data_path <path> -data_name <name> -num_al <number> -ref_panel <path> [-num_cpus <number>] [-num_gpus <number>]"
    exit 1
}

# Parse command line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -data_path) data_path="$2"; shift ;;
        -data_name) data_name="$2"; shift ;;
        -num_al) num_al="$2"; shift ;;
        -ref_panel) ref_panel="$2"; shift ;;
        -num_cpus) num_cpus="$2"; shift ;;
        -num_gpus) num_gpus="$2"; shift ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Check if required parameters are provided
if [ -z "$data_path" ] || [ -z "$data_name" ]  || [ -z "$num_al" ]; then
    usage
fi

# Log the input parameters
echo "Running with the following parameters:"
echo "Data path: $data_path"
echo "Data name: $data_name"
echo "Number of alleles: $num_al"
echo "Reference panel: $ref_panel"
echo "CPUs: $num_cpus"
echo "GPUs: $num_gpus"

# Create output directory
output=${data_path}/output
mkdir -p ${output}

log=$(readlink -f "${output}/${data_name}".log)
echo "Log file: ${log}"

echo "Remove output files if existing:"
if [ -e "${output}"/"${data_name}".imputed.hap ]
then
  echo "- rm ${output}/${data_name}.imputed.hap"
  rm "${output}"/"${data_name}".imputed.hap
fi
if [ -e "${data_path}"/"${data_name}".hap ]
then
  echo "- rm ${data_path}/${data_name}.hap"
  rm "${data_path}"/"${data_name}".hap
fi
if [ -e "${data_path}"/"${data_name}".legend ]
then
  echo "- rm ${data_path}/${data_name}.legend"
  rm "${data_path}"/"${data_name}".legend
fi
if [ -e "${output}"/"${data_name}".imputed.vcf ]
then
  echo "- rm ${output}/${data_name}.imputed.vcf"
  rm "${output}"/"${data_name}".imputed.vcf
fi
if [ -e "${output}"/"${data_name}".imputed.dip ]
then
  echo "- rm ${output}/${data_name}.imputed.dip"
  rm "${output}"/"${data_name}".imputed.dip
fi

echo "Convert diplotype to haplotype: "
echo "- bcftools query -f '[%GT ]\n' ${data_path}/${data_name}.vcf | sed -E 's/[|\/]/ /g' | sed -E 's/\./?/g' | sed 's/ *$//' > ${data_path}/${data_name}.hap"
bcftools query -f '[%GT ]\n' "${data_path}"/"${data_name}".vcf | sed -E 's/[|\/]/ /g' | sed -E 's/\./?/g' | sed 's/ *$//' > "${data_path}"/"${data_name}".hap
echo "- bcftools query -f '[%GT ]\n' "${data_path}"/"${data_name:0:-8}".ori.vcf | sed -E 's/[|\/]/ /g' | sed -E 's/\./?/g' | sed 's/ *$//' > "${data_path}"/"${data_name:0:-8}".ori.hap"
bcftools query -f '[%GT ]\n' "${data_path}"/"${data_name:0:-8}".ori.vcf | sed -E 's/[|\/]/ /g' | sed -E 's/\./?/g' | sed 's/ *$//' > "${data_path}"/"${data_name:0:-8}".ori.hap

echo "Create legend file: "
echo "- bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT' ${data_path}/${data_name}.vcf  > ${data_path}/${data_name}.legend"
bcftools query -f '%CHROM\t%POS\t%ID\t%REF\t%ALT\n' "${data_path}"/"${data_name}".vcf  > "${data_path}"/"${data_name}".legend

if [ -n $ref_panel ]
then
  echo "Find intersection and convert diplotype reference panel to haplotype one: "
  # Define the reference panel file and its gzipped version
  gz="${data_path}/${data_name}.vcf.gz"

  # Check if the gzipped file exists
  if [ ! -f "$gz" ]; then
      # If the gzipped file doesn't exist, compress the original file with bgzip
      bgzip -c "${data_path}/${data_name}.vcf" > "$gz"
      echo "Compressed ${data_path}/${data_name}.vcf to ${gz}."
  else
      # If the gzipped file exists, do nothing
      echo "$gz already exists. Skipping compression."
  fi

  # Index the gzipped file with tabix if the index doesn't already exist
  if [ ! -f "${gz}.tbi" ]; then
      # If the index doesn't exist, create it
      tabix -p vcf "$gz"
      echo "Indexed ${gz} with tabix."
  else
    # If the index exists, do nothing
    echo "${gz}.tbi already exists. Skipping indexing."
  fi

  echo "-  bcftools isec -n=2 -w1 -p ${data_path} ${ref_panel} ${gz} | \
  bcftools query -f '[%GT ]\n' | \
  sed -E 's/[|\/]/ /g' | sed -E 's/\./?/g' | sed 's/ *$//' > ${data_path}/ref_panel.hap"

  bcftools isec -n=2 -w1 -o ${data_path}/isec.vcf ${ref_panel} ${gz}
  bcftools query -f '[%GT ]\n' ${data_path}/isec.vcf | \
  sed -E 's/[|\/]/ /g' | sed -E 's/\./?/g' | sed 's/ *$//' > ${data_path}/ref_panel.hap
  cat ${data_path}/ref_panel.hap ${data_path}/${data_name}.hap > ${data_path}/input.hap
fi

cp ${data_path}/${data_name}.hap  ${data_path}/input.hap

echo "Run imputation with reference panel: ${ref_panel} "
echo """
python3 main.py
--missing_data ${data_path}/${data_name}.hap
--output_data ${output}/${data_name}.imputed.hap
--num_al ${num_al}
--batch_size 100
--hint_rate 0.9
--alpha 100
--iterations 1000
--num_cpus ${num_cpus}
--num_gpus ${num_gpus}
"""
SECONDS=0
python3 main.py \
--missing_data ${data_path}/input.hap \
--output_data ${output}/output.hap \
--num_al ${num_al} \
--batch_size 100 --hint_rate 0.9 --alpha 100 \
--iterations 1000 \
--num_cpus "${num_cpus}" \
--num_gpus "${num_gpus}"

# Count the number of rows (samples) in input
num_samples=$(wc -l < ${data_path}/${data_name}.hap)
# Take the last number of samples from output
tail -n "$num_samples" ${output}/output.hap > ${output}/${data_name}.imputed.hap

"${output}"/"${data_name}".imputed.hap
echo -e "- Sample: ${data_name}\n $(date -d@${SECONDS} -u +%H:%M:%S) seconds elapsed.\n"

echo "Export imputed data:"

echo "- Get header from VCF file: "
echo "  bcftools view -h ${data_path}/${data_name}.vcf | head -n -1 > ${output}/${data_name}.imputed.vcf"
bcftools view -h "${data_path}"/"${data_name}".vcf | head -n -1 > ${output}/${data_name}.imputed.vcf
echo -ne "#CHROM\tPOS\tID\tREF\tALT\t" >> ${output}/${data_name}.imputed.vcf
bcftools query -l ${data_path}/${data_name}.vcf | tr '\n' '\t' >> ${output}/${data_name}.imputed.vcf

echo "- Convert imputed haplotype to diplotype: "
echo "  cat ${output}/${data_name}.imputed.hap | sed -E 's/([0-9?]*)\s([0-9?]*)\s/\1|\2\t/g' > ${output}/${data_name}.imputed.dip"
cat ${output}/${data_name}.imputed.hap | sed -E 's/([0-9?]*)\s([0-9?]*)\s/\1|\2\t/g' > ${output}/${data_name}.imputed.dip

echo "- Add legend and diplotype to file: "
echo "  paste -d'\t' ${data_path}/${data_name}.legend ${output}/${data_name}.imputed.dip >> ${output}/${data_name}.imputed.vcf"
paste -d'\t' ${data_path}/${data_name}.legend ${output}/${data_name}.imputed.dip >> ${output}/${data_name}.imputed.vcf

echo "Execution complete."