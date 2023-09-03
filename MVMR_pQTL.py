from scipy.linalg import lstsq
from numpy import *
import numpy as np
import pandas as pd
import sys
import rpy2.robjects as ro
from rpy2.robjects import r
from rpy2.robjects import pandas2ri

#
# pandas2ri.activate()

"""

getting the causal genes for cases the user cannot supply the LD-matrix.


To run the file, type in the command line python3 run_MVMR.py "/home/user/exposure_outcome.csv"


where: 

"/home/user/exposure_outcome.csv"

sys.argv[1] is the first argument when you run the file and should be the exposure_outcome.csv file 
containing the SNPs to exposure effects and SNPs to outcome effects.

The output after running this file would be a .csv with results from the methods (depending
 on the dimensions, least-squares, generalized method of moments and ratio method)
saved as .csv file in the same directory as the one given for the exposure_outcome.csv file. 

i.e. for this example "/home/user/exposure_outcome_results.csv"

"""

# file_EXEY = sys.argv[1]
# Data = pd.read_csv(file_EXEY, sep=',')
# for (columnName, columnData) in Data.iteritems():
#     if (Data[columnName] == 0).all():
#         Data.drop(columnName, axis=1, inplace=True)
# snps_full = list(Data.loc[:, 'rs_number'].values)
# # Format the list as a single string with elements enclosed in single quotes and separated by commas
# formatted_string = ', '.join([f"'{element}'" for element in snps_full])
#
# # Print the formatted string
# print(formatted_string)
#
# print("#####################", len(snps_full))
# if len(snps_full) != 1:
#     try:
#         cov_data = r("TwoSampleMR::ld_matrix")(snps_full)
#     except Exception as e:
#         error_message = str(e)
#
#     # Extract variant names from the error message
#     error_variants = [line.strip() for line in error_message.split('\n') if line.strip() != '']
#
#     # Create a boolean mask to filter out rows with variant names from the error
#     mask = ~Data['rs_number'].isin(error_variants)
#
#     # Apply the mask to your DataFrame to remove the rows
#     Data = Data[mask]
#     snps = list(Data.loc[:, 'rs_number'].values)
#     cov_data = r("TwoSampleMR::ld_matrix")(snps_full)
#     pandas2ri.activate()
#     pd_from_r_df = ro.conversion.rpy2py(cov_data)
#     cov_data = pd.DataFrame(pd_from_r_df, columns=r('colnames')(cov_data), index=r('colnames')(cov_data), dtype='float')
#     upper_tri = cov_data.where(np.triu(np.ones(cov_data.shape), k=1).astype(bool))
#     to_drop = [column for column in upper_tri.columns if any(upper_tri[column] == 1.0)]
#     cov_data = cov_data.drop(labels=to_drop, axis=1)
#     cov_data = cov_data.drop(labels=to_drop, axis=0)
#
#     cov_EE = cov_data.values
#     print(" det cov data1", np.linalg.det(cov_EE))
#     snps_new = [i.split('_', 1)[0] for i in r('colnames')(cov_data)]
#     # snps_new = [i for i in r('colnames')(cov_data)]
#     cov_data_pruned = pd.DataFrame(cov_EE, columns=snps_new, index=snps_new, dtype='float')
#     snps_Data = [i.split('_', 1)[0] for i in Data.rs_number]
#     Data.rs_number = snps_Data
#     Data = Data[Data['rs_number'].isin(cov_data_pruned.index.values)]
#     Data.reset_index(drop=True, inplace=True)
#
#     snps_final = Data.loc[:, 'rs_number'].values
#     no_snps = len(snps_final)
#     no_genes = len(Data.columns) - 2
#     df = Data.loc[:, ~Data.columns.isin(['rs_number', 'outcome'])]
#     gene_names = column_names = list(df.columns.values)
#     covEY = Data.loc[:, 'outcome'].values
#     covEX = df.values
#     if no_snps == no_genes:
#         if no_snps == 1:
#             b_est = covEY / covEX
#             d1 = {'gene': gene_names,
#                   'Causal Estimate Ratio method': b_est}
#             df1 = pd.DataFrame(data=d1)
#             df1 = df1.set_index('gene')
#         else:
#             b_est = np.linalg.solve(covEX, covEY)
#             d1 = {'gene': gene_names,
#                   'Causal Estimate': b_est}
#             df1 = pd.DataFrame(data=d1)
#             df1 = df1.set_index('gene')
#     else:
#         if no_snps > no_genes:
#             b_est2, res, rnk, sy = lstsq(covEX, covEY)
#             b_est1 = np.linalg.inv(covEX.T @ np.linalg.inv(cov_EE) @ covEX) @ (covEX.T @ np.linalg.inv(cov_EE) @ covEY)
#             d1 = {'gene': gene_names,
#                   'Causal Estimate Least Squares': b_est2, 'Causal Estimate GMM': b_est1}
#             df1 = pd.DataFrame(data=d1)
#             df1 = df1.set_index('gene')
#         else:
#             sys.exit('Error Message : You require at least as many instruments as exposures to run this analysis.')
#
#
# else:
#     no_genes = len(Data.columns) - 2
#     if len(snps_full) >= no_genes:
#         df = Data.loc[:, ~Data.columns.isin(['rs_number', 'outcome'])]
#         gene_names = column_names = list(df.columns.values)
#         covEY = Data.loc[:, 'CAD'].values
#         covEX = df.values
#         b_est = covEY / covEX
#         d1 = {'gene': gene_names,
#               'Causal Estimate Ratio method': b_est[0]}
#         df1 = pd.DataFrame(data=d1)
#         df1 = df1.set_index('gene')
#     else:
#         sys.exit('Error Message : You require at least as many instruments as exposures to run this analysis.')
#
# df1.to_csv(file_EXEY + "_results.csv", sep=",", float_format='%g')



from scipy.linalg import lstsq
from numpy import *
import numpy as np
import pandas as pd
import sys


"""
getting the causal genes for cases the user can supply the LD-matrix.


To run the fil, type in the command line python3 run_MVMR.py "/home/user/exposure_outcome.csv" "/home/user/ld.csv"


where: 

"/home/user/exposure_outcome.csv"

sys.argv[1] is the first argument when you run the file and should be the exposure_outcome.csv file 
containing the SNPs to exposure effects and SNPs to outcome effects.


"/home/user/ld.csv"

sys.argv[2] is the second argument i.e. ld.csv file for the LD matrix.
Please make sure the ordering of the SNPs is same as in the exposure_outcome file


The output after running this file would be a .csv with results from the methods (depending
 on the dimensions, least-squares, generalized method of moments and ratio method)
saved as .csv file in the same directory as the one given for the exposure_outcome.csv file. 

i.e. for this example "/home/user/exposure_outcome_results.csv.csv"

"""

file_EXEY = sys.argv[1]
file_EE = sys.argv[2]
Data = pd.read_csv(file_EXEY, sep=',')
cov_EE = pd.read_csv(file_EE, delimiter=',')

snps = cov_EE.columns.tolist()
cov_EE = cov_EE.values
# cov_data = pd.read_csv(file_EE, delimiter=',')
# print(cov_data.columns)
# print(cov_data.index)
for (columnName, columnData) in Data.iteritems():
    if (Data[columnName] == 0).all():
        Data.drop(columnName, axis=1, inplace=True)
snps_new = [i.split('_', 1)[0] for i in snps]
Data = Data[Data['rs_number'].isin(snps_new)]
Data.reset_index(drop=True, inplace=True)

if len(snps_new) != 1:
    cov_data = pd.DataFrame(cov_EE, columns=snps_new, index=snps_new, dtype='float')
    upper_tri = cov_data.where(np.triu(np.ones(cov_data.shape), k=1).astype(bool))
    to_drop = [column for column in upper_tri.columns if any(upper_tri[column] == 1.0)]
    cov_data = cov_data.drop(labels=to_drop, axis=1)
    cov_data = cov_data.drop(labels=to_drop, axis=0)
    Data = Data[Data.rs_number.isin(cov_data.index.values)]
    Data.reset_index(drop=True, inplace=True)
    if 1E-2 > np.linalg.det(cov_EE) >= 1E-20:
        # losst = ['rs4970834_T_C', 'rs611917_G_A']
        upper_tri = cov_data.where(np.triu(np.ones(cov_data.shape), k=1).astype(bool))
        to_drop = [column for column in upper_tri.columns if
                   any(upper_tri[column] >= float(0.6))]
        cov_data = cov_data.drop(labels=to_drop, axis=1)
        cov_data = cov_data.drop(labels=to_drop, axis=0)
        cov_EE = cov_data.values
        print("det_new1", np.linalg.det(cov_EE))
    elif np.linalg.det(cov_EE) < 1E-20:
        upper_tri = cov_data.where(np.triu(np.ones(cov_data.shape), k=1).astype(bool))
        to_drop = [column for column in upper_tri.columns if
                   any(upper_tri[column] >= float(0.2))]
        cov_data = cov_data.drop(labels=to_drop, axis=1)
        cov_data = cov_data.drop(labels=to_drop, axis=0)
        cov_EE = cov_data.values
        print("det_new2", np.linalg.det(cov_EE))
        print(cov_EE)


    cov_EE = cov_data.values
    print(" det cov data1", np.linalg.det(cov_EE))
    Data = Data[Data.rs_number.isin(cov_data.index.values)]
    Data.reset_index(drop=True, inplace=True)
    print("############ Data", Data)
    snps_new = Data.loc[:, 'rs_number'].values
    no_snps = len(snps_new)
    print("############", snps_new)
    no_genes = len(Data.columns) - 2
    df = Data.loc[:, ~Data.columns.isin(['rs_number', 'outcome'])]
    gene_names = column_names = list(df.columns.values)
    covEY = Data.loc[:, 'outcome'].values
    covEX = df.values
    if no_snps == no_genes:
        if no_snps == 1:
            b_est = covEY / covEX
            d1 = {'gene': gene_names,
                  'Causal Estimate Ratio method': b_est}
            df1 = pd.DataFrame(data=d1)
            df1 = df1.set_index('gene')
        else:
            b_est = np.linalg.solve(covEX, covEY)
            d1 = {'gene': gene_names,
                  'Causal Estimate': b_est}
            df1 = pd.DataFrame(data=d1)
            df1 = df1.set_index('gene')
    else:
        if no_snps > no_genes:
            b_est2, res, rnk, sy = lstsq(covEX, covEY)
            b_est1 = np.linalg.inv(covEX.T @ np.linalg.inv(cov_EE) @ covEX) @ (covEX.T @ np.linalg.inv(cov_EE) @ covEY)
            d1 = {'gene': gene_names,
                  'Causal Estimate Least Squares': b_est2, 'Causal Estimate GMM': b_est1}
            df1 = pd.DataFrame(data=d1)
            df1 = df1.set_index('gene')
        else:
            sys.exit('Error Message : You require at least as many instruments as exposures to run this analysis.')

else:
    no_genes = len(Data.columns) - 2
    if len(snps) >= no_genes:
        df = Data.loc[:, ~Data.columns.isin(['SNPs', 'outcome'])]
        gene_names = column_names = list(df.columns.values)
        covEY = Data.loc[:, 'outcome'].values
        covEX = df.values
        b_est = covEY / covEX
        d1 = {'gene': gene_names,
              'Causal Estimate Ratio method': b_est[0]}
        df1 = pd.DataFrame(data=d1)
        df1 = df1.set_index('gene')
    else:
        sys.exit('Error Message : You require at least as many instruments as exposures to run this analysis.')
df1.to_csv(file_EXEY + "_results.csv", sep=",", float_format='%g')
