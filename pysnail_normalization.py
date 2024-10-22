import pandas
import os

def pysnail(expression_data_path, group_data_path, output_path):
    """
    Apply quantile normalization to expression data.

    This function reads expression data and group information from the provided
    file paths, performs quantile smooth normalization on the expression data, and
    saves the normalized data to the specified output path.

    Parameters:
    expression_data_path (str): Path to the input expression data file in TSV format.
    group_data_path (str): Path to the input group data file in TSV format.
    output_path (str): Path to save the normalized expression data file in TSV format.

    Returns:
    pandas.DataFrame: DataFrame containing the normalized expression data.
    """
    xprs = os.path.realpath(expression_data_path)
    groups = os.path.realpath(group_data_path)
    output_file = os.path.realpath(output_path)

    dataset = Dataset(xprs, groups, **{'index_col': 0, 'sep':'\t'})
    xprs_norm, qstat = qsmooth(dataset, aggregation='auto', cutoff=0.15)

    xprs_norm.to_csv(output_file, sep='\t')
    print(f"Normalized expression data saved to {output_file}")
    
    return xprs_norm

if __name__ == "__main__":
