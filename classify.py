import os
import pickle
import argparse
import pandas as pd
import extract_features


def classify_data(model, test_folder, refseq_database):
    # Extract features from the training folder.
    if not os.path.isfile(os.path.join(test_folder, 'extracted_features.csv')):
        extract_features.main(sequencepath=test_folder,
                              report=True,
                              refseq_database=refseq_database)
    test_df = pd.read_csv(os.path.join(test_folder, 'extracted_features.csv'))
    dataframe = pd.get_dummies(test_df, columns=['Genus'], dummy_na=True)
    features = list(dataframe.columns[1:len(dataframe.columns)])
    x = dataframe[features]
    result = model.predict(x)
    # result = tree.predict_proba(x)
    for i in range(len(result)):
        if result[i] == 0:
            output = 'Fail'
        elif result[i] == 1:
            output = 'Pass'
        elif result[i] == 2:
            output = 'Reference'
        print(test_df['SampleName'][i] + ',' + output)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--test_folder',
                        type=str,
                        required=True,
                        help='Path to folder containing FASTA files you want to test.')
    parser.add_argument('-m', '--model',
                        type=str,
                        required=True,
                        help='Path to pickle object of model you want to use for classification.')
    parser.add_argument('-d', '--refseq_database',
                        type=str,
                        required=True,
                        help='Path to custom refseq database.')
    args = parser.parse_args()
    random_forest_model = pickle.load(open(args.model, 'rb'))
    print(sorted(zip(map(lambda x: round(x, 4), random_forest_model.feature_importances_))))
    classify_data(random_forest_model, args.test_folder, args.refseq_database)
