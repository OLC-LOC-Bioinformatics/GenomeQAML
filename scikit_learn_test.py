import os
import click
import sklearn
import subprocess
import numpy as np
import pandas as pd
import extract_features
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split, cross_val_score


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
    dataframe = pd.get_dummies(dataframe, columns=['Genus'], dummy_na=True)
    features = list(dataframe.columns[1:len(dataframe.columns)])
    features.remove('PassFail')
    X = dataframe[features]
    y = dataframe['PassFail']
    X_train, X_test, y_train, y_test = train_test_split(X, y)
    dt = DecisionTreeClassifier()
    # dt = RandomForestClassifier()
    scores = cross_val_score(dt, X, y, cv=5)
    print(np.mean(scores))
    dt = dt.fit(X_train, y_train)
    with open('dt.dot', 'w') as f:
        sklearn.tree.export_graphviz(dt, out_file=f, feature_names=features)
    command = ["dot", "-Tpng", "dt.dot", "-o", "dt.png"]
    subprocess.check_call(command)
    return dt


def predict_results(fasta_dir, tree):
    test_df = pd.read_csv(os.path.join(fasta_dir, 'extracted_features.csv'))
    dataframe = pd.get_dummies(test_df, columns=['Genus'], dummy_na=True)
    features = list(dataframe.columns[1:len(dataframe.columns)])
    # features = list(test_df.columns[1:23])
    x = dataframe[features]
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
@click.option('-d', '--refseq_database',
              type=click.Path(exists=True),
              required=True,
              help='Path to reduced refseq database sketch.')
def cli(pass_folder, fail_folder, test_folder, refseq_database):
    extract_features.main(sequencepath=fail_folder,
                          refseq_database=refseq_database,
                          report=True)
    extract_features.main(sequencepath=pass_folder,
                          refseq_database=refseq_database,
                          report=True)
    df = combine_csv_files(fail_folder=fail_folder, pass_folder=pass_folder)
    # extract_features.main(sequencepath=test_folder,
    #                       refseq_database=refseq_database,
    #                       report=True)
    dt = fit_model(df)
    # predict_results(test_folder, dt)


if __name__ == '__main__':
    cli()
