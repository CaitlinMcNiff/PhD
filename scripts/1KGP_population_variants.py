import os
import gzip
import matplotlib.pyplot as plt

# Initialize counters for each population
populations = ['EAS', 'AMR', 'AFR', 'EUR', 'SAS']
pop_count = {pop: 0 for pop in populations}
unique_pop_count = {pop: 0 for pop in populations}

# Initialise file for storing variants for each populations
variant_files = {pop: open(f'{pop}_variants.vcf', 'w') for pop in populations}

# Function to process a single file
def process_file(filepath):
    with gzip.open(filepath, 'rt') as f:
        for line in f:
            if line.startswith('#'):
                continue

            # Split the line by tabs and defined the desired columns
            cols = line.strip().split('\t')
            chr = cols[0]
            start_pos = int(cols[1])
            end_pos = start_pos + 1
            rsid = cols[2] # rsID
            ra = cols[3] # reference allele
            aa = cols[4] # alternative allele
            info_field = cols[7] # allele frequency


            # Extract allele frequency data for each population
            eas_af = extract_af(info_field, 'EAS_AF=')
            amr_af = extract_af(info_field, 'AMR_AF=')
            afr_af = extract_af(info_field, 'AFR_AF=')
            eur_af = extract_af(info_field, 'EUR_AF=')
            sas_af = extract_af(info_field, 'SAS_AF=')

            # Create a list of AF values
            af_values = {
                'EAS': eas_af,
                'AMR': amr_af,
                'AFR': afr_af,
                'EUR': eur_af,
                'SAS': sas_af
            }

            # Count each variant in the populations it appears in
            populations_present = [pop for pop, af in af_values.items() if af > 0]
            for pop in populations_present:
                pop_count[pop] += 1 
                # write variant information to the population-specific file
                variant_files[pop].write(f"{chr}\t{start_pos}\t{end_pos}\t{rsid}\t{ra}\t{aa}\t{af_values[pop]}\n")

            # Count variants found in only one population
            if len(populations_present) == 1:
                unique_pop_count[populations_present[0]] += 1

# Helper function to extract allele frequency from the info field
def extract_af(info_field, af_key):
    af_info = [x for x in info_field.split(';') if x.startswith(af_key)]
    if af_info:
        # Handle cases where there are multiple AF values (e.g., '0.002,0.001')
        af_values = af_info[0].split('=')[1].split(',')
        af_values = [float(af) for af in af_values if af != '.']  # Convert to float, ignoring missing values
        if af_values:
            return max(af_values)  # Use the max AF value
    return 0.0


# Function to process all files in the folder
def process_all_files_in_folder(folder_path):
    for filename in os.listdir(folder_path):
        if filename.endswith('.vcf.gz'):
            filepath = os.path.join(folder_path, filename)
            print(f'Processing file: {filename}')
            process_file(filepath)

# Function to save data to a text file
def save_data_to_file(data, filename):
    with open(filename, 'w') as f:
        for pop, count in data.items():
            f.write(f'{pop}: {count}\n')

# Main execution
if __name__ == '__main__':
    folder_path = '/exports/cmvm/eddie/sbms/groups/young-lab/caitlin/phd/1KGP_hg38/'
    process_all_files_in_folder(folder_path)

    # Save pop_count and unique_pop_count to text files
    save_data_to_file(pop_count, 'pop_count.txt')
    save_data_to_file(unique_pop_count, 'unique_pop_count.txt')

    # close the variant files for each population
    for file in variant_files.values():
        file.close()