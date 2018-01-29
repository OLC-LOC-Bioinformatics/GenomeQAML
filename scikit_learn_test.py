import os
import click
import pickle
import numpy as np
import pandas as pd
import extract_features
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import cross_val_score


def combine_csv_files(pass_folder, fail_folder, ref_folder):
    df_fail = pd.read_csv(os.path.join(fail_folder, 'extracted_features.csv'))
    df_pass = pd.read_csv(os.path.join(pass_folder, 'extracted_features.csv'))
    df_ref = pd.read_csv(os.path.join(ref_folder, 'extracted_features.csv'))
    # Add a column of 0s to fail, 1s to pass.
    fail = [0] * len(df_fail)
    not_fail = [1] * len(df_pass)
    ref = [2] * len(df_ref)
    df_fail['PassFail'] = fail
    df_pass['PassFail'] = not_fail
    df_ref['PassFail'] = ref
    frames = [df_fail, df_pass, df_ref]
    result = pd.concat(frames)
    return result


def fit_model(dataframe):
    # Use pandas to get one-hot-encoding of categorical featues (Genus).
    dataframe = pd.get_dummies(dataframe, columns=['Genus'], dummy_na=True)
    features = list(dataframe.columns[1:len(dataframe.columns)])
    features.remove('PassFail')  # Make sure PassFail isn't a feature.
    X = dataframe[features]
    y = dataframe['PassFail']
    dt = RandomForestClassifier(max_depth=10, max_leaf_nodes=20)
    # dt = RandomForestClassifier()
    scores = cross_val_score(dt, X, y, cv=10)
    print(np.mean(scores))
    dt = dt.fit(X, y)
    return dt


def predict_results(fasta_dir, tree, training_dataframe):
    test_df = pd.read_csv(os.path.join(fasta_dir, 'extracted_features.csv'))
    dataframe = pd.get_dummies(test_df, columns=['Genus'], dummy_na=True)
    # Remove any genera from the test dataframe that weren't part of our training set.
    training_dataframe = pd.get_dummies(training_dataframe, columns=['Genus'], dummy_na=True)
    for column in dataframe:
        if column not in training_dataframe:
            dataframe.drop(column, axis=1, inplace=True)
    # Add any genera that weren't in our test set but were in the training set.
    for column in training_dataframe:
        if 'Genus' in column and column not in dataframe:
            not_present = [0] * len(dataframe)
            dataframe[column] = not_present
    # Then, add any features that were part of training data but not part of test data
    features = list(dataframe.columns[1:len(dataframe.columns)])
    x = dataframe[features]
    result = tree.predict(x)
    # result = tree.predict_proba(x)
    for i in range(len(result)):
        if result[i] == 0:
            output = 'Fail'
        elif result[i] == 1:
            output = 'Pass'
        elif result[i] == 2:
            output = 'Reference'
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
@click.option('-r', '--ref_folder',
              type=click.Path(exists=True),
              required=True,
              help='Path to folder with training FASTAs that have excellent quality metrics.')
@click.option('-t', '--test_folder',
              type=click.Path(exists=True),
              required=True,
              help='Path to folder with FASTAs you want to test.')
@click.option('-d', '--refseq_database',
              type=click.Path(exists=True),
              required=True,
              help='Path to reduced refseq database sketch.')
def cli(pass_folder, fail_folder, test_folder, refseq_database, ref_folder):
    # Extract features for pass data, fail data, and reference data if it hasn't already been done.
    if not os.path.isfile(os.path.join(fail_folder, 'extracted_features.csv')):
        extract_features.main(sequencepath=fail_folder,
                              refseq_database=refseq_database,
                              report=True)
    if not os.path.isfile(os.path.join(pass_folder, 'extracted_features.csv')):
        extract_features.main(sequencepath=pass_folder,
                              refseq_database=refseq_database,
                              report=True)
    if not os.path.isfile(os.path.join(ref_folder, 'extracted_features.csv')):
        extract_features.main(sequencepath=ref_folder,
                              refseq_database=refseq_database,
                              report=True)

    # Combine the dataframes for training data so that we can fit our decision tree.
    df = combine_csv_files(fail_folder=fail_folder, pass_folder=pass_folder, ref_folder=ref_folder)
    dt = fit_model(df)

    # Extract features for our test set if it hasn't already been done and attempt to predict results.
    if not os.path.isfile(os.path.join(test_folder, 'extracted_features.csv')):
        extract_features.main(sequencepath=test_folder,
                              refseq_database=refseq_database,
                              report=True)
    predict_results(test_folder, dt, df)  # TODO: Add check that FASTA folder actuall has stuff in it.
    pickle.dump(dt, open('model.p', 'wb'))

if __name__ == '__main__':
    cli()
