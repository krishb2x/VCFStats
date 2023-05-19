import gzip
import os
import sys

vcf_dir = sys.argv[1]

def count_variant_metrics(vcf_file):
    total_variant_count = 0
    passed_variant_count = 0
    snp_count = 0
    insertion_count = 0
    deletion_count = 0
    complex_indel_count = 0
    mixed_count = 0
    het_count = 0
    hom_var_count = 0
    singleton_count = 0
    called_genotype_count = 0

    with gzip.open(vcf_file, 'rt') as file:
        for line in file:
            if line.startswith('#'):
                continue

            fields = line.split('\t')
            alt_alleles = fields[4].split(',')

            total_variant_count += 1

            if 'PASS' in fields[6]:
                passed_variant_count += 1

            if len(alt_alleles) == 1 and len(alt_alleles[0]) == 1:
                snp_count += 1
            elif len(alt_alleles) == 1 and len(alt_alleles[0]) > 1:
                insertion_count += 1
            elif len(alt_alleles) > 1 and any(len(alt) > 1 for alt in alt_alleles):
                complex_indel_count += 1
            elif len(alt_alleles) == 1 and alt_alleles[0] == '*':
                deletion_count += 1
            else:
                mixed_count += 1

            genotype_info = fields[9].split(':')
            gt_field = genotype_info[0]

            if gt_field == '0/1':
                het_count += 1
            elif gt_field == '1/1':
                hom_var_count += 1

            if gt_field != './.':
                called_genotype_count += 1

            if len(alt_alleles) == 1 and alt_alleles[0] != '*':
                if gt_field == '0/1' or gt_field == '1/1':
                    singleton_count += 1

    return total_variant_count, passed_variant_count, snp_count, insertion_count, deletion_count, complex_indel_count, mixed_count, het_count, hom_var_count, singleton_count, called_genotype_count


subdir_list = os.listdir(vcf_dir)

for subdir in subdir_list:
    subdir_path = os.path.join(vcf_dir, subdir)
    if os.path.isdir(subdir_path):
        file_list = os.listdir(subdir_path)

        for file in file_list:
            if file.endswith('.hard-filtered.vcf.gz'):
                vcf_path = os.path.join(subdir_path, file)
                sample_name = file.split('.')[0]
                output_file = f'{sample_name}.variant_metrics.txt'

                total_variant, passed_variant, snp, insertion, deletion, complex_indel, mixed, het, hom_var, singleton, called_genotypes = count_variant_metrics(vcf_path)

                # Save the output in a text file
                output_path = os.path.join(subdir_path, output_file)
                with open(output_path, 'w') as output_file_obj:
                    output_file_obj.write(f'Sample Name: {sample_name}\n')
                    output_file_obj.write(f'Number of total variant loci: {total_variant}\n')
                    output_file_obj.write(f'Number of passed variant loci: {passed_variant}\n')
                    output_file_obj.write(f'Number of SNP loci: {snp}\n')
                    output_file_obj.write(f'Number of insertions: {insertion}\n')
                    output_file_obj.write(f'Number of deletions: {deletion}\n')
                    output_file_obj.write(f'Number of complex indels: {complex_indel}\n')
                    output_file_obj.write(f'Number of mixed loci: {mixed}\n')
                    output_file_obj.write(f'Number of het loci: {het}\n')
                    output_file_obj.write(f'Number of hom var loci: {hom_var}\n')
                    output_file_obj.write(f'Number of singletons: {singleton}\n')
                    output_file_obj.write(f'Number of called genotypes: {called_genotypes}\n')

                print(f'Output saved in vcf_path + {output_path}')
