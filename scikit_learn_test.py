import os
import click
import sklearn
import subprocess
import numpy as np
import pandas as pd
import extract_features
from sklearn.tree import DecisionTreeClassifier
from sklearn.model_selection import train_test_split, cross_val_score


def generate_pass_csv(pass_folder):  # This method is now somewhat pointless.
    extract_features.main(pass_folder, report=True)


def generate_fail_csv(fail_folder):
    extract_features.main(fail_folder, report=True)


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
    features = list(dataframe.columns[1:22])
    X = dataframe[features]
    y = dataframe['PassFail']
    X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=0)
    dt = DecisionTreeClassifier()
    scores = cross_val_score(dt, X, y, cv=5)
    print(np.mean(scores))
    dt = dt.fit(X_train, y_train)
    with open('dt.dot', 'w') as f:
        sklearn.tree.export_graphviz(dt, out_file=f, feature_names=features)
    command = ["dot", "-Tpng", "dt.dot", "-o", "dt.png"]
    subprocess.check_call(command)
    return dt


def predict_results(fasta_dir, tree):
    extract_features.main(fasta_dir, report=True)
    test_df = pd.read_csv(os.path.join(fasta_dir, 'extracted_features.csv'))
    print(test_df)
    features = list(test_df.columns[1:22])
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
