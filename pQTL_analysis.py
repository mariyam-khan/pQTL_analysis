import numpy as np
import pandas as pd
import os
from collections import Counter

data_generef = pd.read_csv('~/pQTL_MVMR/pQTL_MVMR/gene_ref.tsv', delimiter='\t')
# Loop through the rows of the 'Gene stable ID' column
DATA2 = pd.read_csv('/Home/ii/mariyamk/pQTL_MVMR/pQTL_MVMR/SNP_ref.tsv', delimiter='\t')
DATA3 = pd.read_csv(
    '/Home/ii/mariyamk/pQTL_MVMR/pQTL_MVMR/plasma_STARNET_cases_matched_cis-pQTLs_500Kb.tsv',
    delimiter='\t')
for i in range(len(data_generef)):
    # Access the first row and first column element
    gene = data_generef.iloc[i, 0]
    chr = data_generef.iloc[i, 1]
    tss = data_generef.iloc[i, 2]

    # Calculate the interval
    interval_size = 500000  # 500KB

    data_generef['interval_start'] = data_generef['tss'] - interval_size
    data_generef['interval_end'] = data_generef['tss'] + interval_size

    target_gene_chromosome = chr
    target_gene_start = data_generef[data_generef['Gene stable ID'] == gene]['interval_start'].values[0]
    target_gene_end = data_generef[data_generef['Gene stable ID'] == gene]['interval_end'].values[0]

    # Now, filter the rs_numbers in DATA2 that match the chromosome and fall within the interval
    filtered_rs_numbers = DATA2[
        (DATA2['chromosome'] == target_gene_chromosome) & (DATA2['position'] >= target_gene_start) & (
                DATA2['position'] <= target_gene_end)]['rs_number']

  
    list = ['AOR', 'Blood', 'LIV', 'MAM', 'SF', 'SKLM', 'VAF']
    # Check if filtered_rs_numbers is not empty
    if not filtered_rs_numbers.empty:
        print("################################## filtered_rs_numbers", filtered_rs_numbers)
        # Define the list of suffixes
        suffix_list = ['AOR', 'Blood', 'LIV', 'MAM', 'SF', 'SKLM', 'VAF']

        # Define the directory where your files are located
        directory_path = '/Home/ii/mariyamk/pQTL_MVMR/pQTL_MVMR'  # Replace with the actual path

        # Create an empty DataFrame to store the collected data
        collected_data = pd.DataFrame()

        # Loop through the files with specified suffixes
        for suffix in suffix_list:
            # Construct the file name based on the suffix
            file_name = f'{suffix}_STARNET_cases_matched_QTLs_cis-eQTLs_500Kb.tsv'  # Replace with your file naming convention

            # Construct the full file path
            file_path = os.path.join(directory_path, file_name)

            # Check if the file exists
            if os.path.exists(file_path):
                # Load the file into a DataFrame
                df = pd.read_csv(file_path, delimiter='\t')  # You may need to adjust the read_csv parameters

                # Filter the data based on filtered_rs_numbers
                filtered_data = df[df['rs_number'].isin(filtered_rs_numbers)]

                num_rows = len(filtered_data)
                print("num_rows", num_rows)
                while num_rows > 500:
                    interval_size -= 1000  # Decrease interval size
                    print("Length greater than 500. Reducing interval size to:", interval_size)
                    data_generef['interval_start'] = data_generef['tss'] - interval_size
                    data_generef['interval_end'] = data_generef['tss'] + interval_size

                    target_gene_start = data_generef[data_generef['Gene stable ID'] == gene]['interval_start'].values[0]
                    target_gene_end = data_generef[data_generef['Gene stable ID'] == gene]['interval_end'].values[0]

                    # Now, filter the rs_numbers in DATA2 that match the chromosome and fall within the interval
                    filtered_rs_numbers = DATA2[
                        (DATA2['chromosome'] == target_gene_chromosome) & (DATA2['position'] >= target_gene_start) & (
                                DATA2['position'] <= target_gene_end)]['rs_number']

                    filtered_data = df[df['rs_number'].isin(filtered_rs_numbers)]
                    filtered_data.reset_index(drop=True)
                    num_rows = len(filtered_data)
                    print("num_rows_new", num_rows)
                # Check if filtered_data is not empty
                if not filtered_data.empty:
                    # Add a new column for gene name + suffix
                    filtered_data['exposure'] = gene + '_' + suffix
                    # Append the filtered data to the collected_data DataFrame
                    collected_data = collected_data.append(filtered_data, ignore_index=True)
        # Drop and reset the index
        collected_data = collected_data.reset_index(drop=True)
        collected_data = collected_data.rename(columns={'beta': 'beta.exposure'})
        # Now, collected_data contains the data from all the files that match filtered_rs_numbers
        # Get unique rs_number values from collected_data
        unique_rs_numbers = collected_data['rs_number'].unique()

        # Create a new DataFrame to store the results
        result_df = pd.DataFrame({'rs_number': unique_rs_numbers})
        # Merge result_df with DATA3 to get the beta values
        result_df = result_df.merge(DATA3[['rs_number', 'beta', 'p-value']], on='rs_number', how='left')

        # Rename the 'beta' column to 'beta.outcome'
        result_df = result_df.rename(columns={'beta': 'beta.outcome'})
       
        Data_out = result_df
        Data_out = Data_out.sort_values(by='p-value', ascending=True)
        Data_out.reset_index(drop=True, inplace=True)
        Data_exp = collected_data

        unique_genes = Counter(Data_exp['exposure']).keys()
        no_genes = len(unique_genes)
        most_common = Counter(Data_exp['exposure']).most_common(no_genes)
        genes = []
        for j in range(no_genes):
            genes.append(most_common[j][0])
        selected_columns2 = Data_out['beta.outcome']
        new_df2 = selected_columns2.copy()
        column_names = genes
        row_names = unique_rs_numbers
        matrix = np.zeros((len(row_names), len(column_names)))
        R_data = pd.DataFrame(matrix, columns=column_names, index=row_names, dtype='object')
        s_data = pd.DataFrame(np.zeros(len(row_names)), columns=['outcome'], index=row_names, dtype='object')

        for n in range(no_genes):
            name = genes[n]
            exp = (most_common[n][0])
            Data = Data_exp[Data_exp['exposure'] == exp]
            Data.reset_index(drop=True, inplace=True)
            o = Data['rs_number'].values
            for j1 in range(int(most_common[n][1])):
                value = o[j1]
                R_data.at[value, name] = Data.loc[j1, 'beta.exposure']
                idx = Data_out[Data_out['rs_number'] == value].index.values
                s_data.loc[value, 'outcome'] = new_df2[idx[0]]

        R_data.replace(np.nan, 0)
        s_data.replace(np.nan, 0)
        df_R = pd.DataFrame(data=R_data)
        df_s = pd.DataFrame(data=s_data)
        RandS = pd.concat([df_R, df_s], axis=1)
        RandS.index.name = 'rs_number'
        RandS.to_csv("/Home/ii/mariyamk/pQTL_MVMR/pQTL_MVMR/data_" + gene + ".csv", sep=",")
        print("########### saved ")
    else:
        print("filtered_rs_numbers is empty.")
