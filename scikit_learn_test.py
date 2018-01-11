import click
import os
import subprocess
import pandas as pd
import sklearn
from sklearn.tree import DecisionTreeClassifier
import extract_features


def generate_pass_csv(pass_folder):
    gc_dict, longest_contig_dict, genome_length_dict, num_contigs_dict, n50_dict, \
        n75_dict, l50_dict, l75_dict = extract_features.main(pass_folder)
    extract_features.report(gc_dict, longest_contig_dict, genome_length_dict, num_contigs_dict, n50_dict,
                            n75_dict, l50_dict, l75_dict, pass_folder)


def generate_fail_csv(fail_folder):
    gc_dict, longest_contig_dict, genome_length_dict, num_contigs_dict, n50_dict, \
        n75_dict, l50_dict, l75_dict = extract_features.main(fail_folder)
    extract_features.report(gc_dict, longest_contig_dict, genome_length_dict, num_contigs_dict, n50_dict,
                            n75_dict, l50_dict, l75_dict, fail_folder)


def combine_csv_files(pass_folder, fail_folder):
    df_fail = pd.read_csv(os.path.join(fail_folder, 'extracted_features.csv'))
    df_pass = pd.read_csv(os.path.join(pass_folder, 'extracted_features.csv'))
    # Add a column of 0s to fail, 1s to pass.
    fail = [0] * len(df_fail)
    not_fail = [1] * len(df_pass)
    df_fail['PassFail'] = fail
    df_pass['PassFail'] = not_fail
    frames = [df_fail, df_pass]
    result = pd.concat(frames)
    return result


def fit_model(dataframe):
    features = list(dataframe.columns[1:9])
    y = dataframe['PassFail']
    x = dataframe[features]
    dt = DecisionTreeClassifier()
    dt = dt.fit(x, y)
    with open('dt.dot', 'w') as f:
        sklearn.tree.export_graphviz(dt, out_file=f, feature_names=features)
    command = ["dot", "-Tpng", "dt.dot", "-o", "dt.png"]
    subprocess.check_call(command)
    return dt


def predict_results(fasta_dir, tree):
    gc_dict, longest_contig_dict, genome_length_dict, num_contigs_dict, n50_dict, \
        n75_dict, l50_dict, l75_dict = extract_features.main(fasta_dir)
    extract_features.report(gc_dict, longest_contig_dict, genome_length_dict, num_contigs_dict, n50_dict,
                            n75_dict, l50_dict, l75_dict, fasta_dir)
    test_df = pd.read_csv(os.path.join(fasta_dir, 'extracted_features.csv'))
    features = list(test_df.columns[1:9])
    x = test_df[features]
    result = tree.predict(x)
    for i in range(len(result)):
        if result[i] == 0:
            output = 'Fail'
        elif result[i] == 1:
            output = 'Pass'
        print(test_df['SampleName'][i] + ',' + output)


@click.command()
@click.option('-p', '--pass_folder',
              type=click.Path(exists=True),
              required=True,
              help='Path to folder with training FASTAs that pass quality metrics.')
@click.option('-f', '--fail_folder',
              type=click.Path(exists=True),
              required=True,
              help='Path to folder with training FASTAs that fail quality metrics.')
@click.option('-t', '--test_folder',
              type=click.Path(exists=True),
              required=True,
              help='Path to folder with FASTAs you want to test.')
def cli(pass_folder, fail_folder, test_folder):
    generate_fail_csv(fail_folder)
    generate_pass_csv(pass_folder)
    df = combine_csv_files(fail_folder=fail_folder, pass_folder=pass_folder)
    dt = fit_model(df)
    predict_results(test_folder, dt)


if __name__ == '__main__':
    cli()