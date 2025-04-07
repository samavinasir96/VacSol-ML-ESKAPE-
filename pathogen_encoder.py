import pandas as pd

def encode_pathogens(request, progress_callback):

    progress_callback("Encoding Selected ESKAPE Pathogen")

    pathogen_encoding = {
        'Staphylococcus aureus': 5,
        'Enterococcus faecium': 2,
        'Pseudomonas aeruginosa': 4,
        'Acinetobacter baumannii': 0,
        'Klebsiella pneumoniae': 3,
        'Enterobacter': 1
    }

    final = pd.read_csv('finalfeatures_updated.csv')

    Protein_ID = final['Protein_ID']
    finalDF = final.drop(labels=["Protein_ID"], axis=1)

    if request.method == 'POST':
        selected_pathogen = request.POST.get('pathogen')
        
        encoded_pathogen = pathogen_encoding.get(selected_pathogen, -1)  # Default to -1 if not found
        
        if encoded_pathogen != -1:

            finalDF["Organism.1"] = encoded_pathogen

            finalDF.to_csv('finalfeatures_pathogen.csv', index=False)
        else:
            print("Error: Selected pathogen not found in the encoding dictionary.")
        