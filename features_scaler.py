import joblib
import pandas as pd


def scale_features(progress_callback):

    progress_callback("Scaling Features")

    final = pd.read_csv('finalfeatures_updated.csv')

    Protein_ID = final['Protein_ID']
    finalDF = final.drop(labels=["Protein_ID"], axis=1)

    scaler_path = 'eskapeml/scaler.pkl'

    dataset = pd.read_csv('finalfeatures_pathogen.csv')
    

    with open(scaler_path, "rb") as scaler:
        scale_dataset = joblib.load(scaler)    
        rescaled=scale_dataset.transform(dataset)
    
    rescaleData = pd.DataFrame(rescaled, columns=finalDF.columns)
    rescaleData.insert(0, 'Protein_ID', Protein_ID)
    rescaleData.to_csv('final_dataset.csv', index=False)

