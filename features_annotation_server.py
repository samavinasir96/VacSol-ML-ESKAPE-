import subprocess
import time
import zipfile
import iedb
import joblib
import pandas as pd
import os
import sys    
from django.http import HttpRequest, HttpResponse, HttpResponseRedirect, JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.shortcuts import redirect, render
import pickle
import pandas as pd
import os
import sys
import logging
import zipfile
from Bio import SeqIO
import biolib


@csrf_exempt
def home(request):
    return render (request, 'index.html')

def get_latest_logs(request):
    script_dir = os.path.dirname(os.path.abspath(__file__))
    log_file_path = os.path.join(script_dir, "vacsol_ml.log")
    with open(log_file_path, 'r') as log_file:
        latest_logs = log_file.readlines()
    return JsonResponse({'logs': latest_logs})

def progress_callback(message, progress):
    logger = logging.getLogger(__name__)
    log_message = f"{message}"
    logger.info(log_message, extra={'progress': progress})

### ESKAPE ###

def calculate_features(request, file_path, progress_callback):
    total_steps = 8
    completed_steps = 0
    # Use os.path.abspath to get the absolute file path
    script_dir = os.path.dirname(os.path.abspath(__file__))
    file_path = os.path.join(script_dir, "sequences.fasta")
    file_path = os.path.abspath(file_path)
    MERGED_DATAFRAMES = []
    protein_ids = []
    mhci_prob_scores = []
    mhci_ranks = []
    mhcii_ranks = []
    mhcii_scores = []
    b_cells_prob_scores = []
    surface_probs = []
    antigenicities = []
    header_list = []

    progress_callback(f"Step 1: Sequence(s) uploaded successfully", 1)
    progress_callback("Sending POST request to Physiochemical parameters analysis tool", 2)

    path_to_virtualenv = '/home/mgbisewu/virtualenv/eskapeml/3.9'
    
    ifeature_path = path_to_virtualenv + '/lib/python3.9/site-packages/iFeature/iFeature.py'
    
    python_path = path_to_virtualenv + '/bin/python3.9'
    
    command1 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type AAC --out biolib_results/ifeature1.tsv'
    command2 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type GAAC --out biolib_results/ifeature2.tsv'
    command3 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type Moran --out biolib_results/ifeature3.tsv'
    command4 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type Geary --out biolib_results/ifeature4.tsv'
    command5 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type NMBroto --out biolib_results/ifeature5.tsv'
    command6 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type CTDC --out biolib_results/ifeature6.tsv'
    command7 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type CTDD --out biolib_results/ifeature7.tsv'
    command8 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type CTDT --out biolib_results/ifeature8.tsv'
    command9 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type CTriad --out biolib_results/ifeature9.tsv'
    command10 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type KSCTriad --out biolib_results/ifeature10.tsv'
    command11 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type SOCNumber --out biolib_results/ifeature11.tsv'
    command12 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type QSOrder --out biolib_results/ifeature12.tsv'
    command13 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type PAAC --out biolib_results/ifeature13.tsv'
    command14 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type APAAC --out biolib_results/ifeature14.tsv'
    command15 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type CKSAAP --out biolib_results/ifeature15.tsv'
    command16 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type DPC --out biolib_results/ifeature16.tsv'
    command17 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type DDE --out biolib_results/ifeature17.tsv'
    command18 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type TPC --out biolib_results/ifeature18.tsv'
    command19 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type CKSAAGP --out biolib_results/ifeature19.tsv'
    command20 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type GDPC --out biolib_results/ifeature20.tsv'
    command21 = f'{python_path} {ifeature_path} --file eskapeml/sequences.fasta --type GTPC --out biolib_results/ifeature21.tsv'

    
    # Run the command
    subprocess.run(command1, shell=True)
    subprocess.run(command2, shell=True)
    subprocess.run(command3, shell=True)
    subprocess.run(command4, shell=True)
    subprocess.run(command5, shell=True)
    subprocess.run(command6, shell=True)
    subprocess.run(command7, shell=True)
    subprocess.run(command8, shell=True)
    subprocess.run(command9, shell=True)
    subprocess.run(command10, shell=True)
    subprocess.run(command11, shell=True)
    subprocess.run(command12, shell=True)
    subprocess.run(command13, shell=True)
    subprocess.run(command14, shell=True)
    subprocess.run(command15, shell=True)
    subprocess.run(command16, shell=True)
    subprocess.run(command17, shell=True)
    subprocess.run(command18, shell=True)
    subprocess.run(command19, shell=True)
    subprocess.run(command20, shell=True)
    subprocess.run(command21, shell=True)
    
    completed_steps += 1
    progress=int((completed_steps / total_steps) * 100)
    progress_callback("Step 2: Physicochemical analysis completed", 3)
    
    df1 = pd.read_table ('biolib_results/ifeature1.tsv')
    df2 = pd.read_table ('biolib_results/ifeature2.tsv')
    df3 = pd.read_table ('biolib_results/ifeature3.tsv')
    df4 = pd.read_table('biolib_results/ifeature4.tsv')
    df5 = pd.read_table('biolib_results/ifeature5.tsv')
    df6 = pd.read_table ('biolib_results/ifeature6.tsv')
    df7 = pd.read_table ('biolib_results/ifeature7.tsv')
    df8 = pd.read_table ('biolib_results/ifeature8.tsv')
    df9 = pd.read_table ('biolib_results/ifeature9.tsv')
    df10 = pd.read_table('biolib_results/ifeature10.tsv')
    df11 = pd.read_table('biolib_results/ifeature11.tsv')
    df12 = pd.read_table ('biolib_results/ifeature12.tsv')
    df13 = pd.read_table ('biolib_results/ifeature13.tsv')
    df14 = pd.read_table ('biolib_results/ifeature14.tsv')
    df15 = pd.read_table('biolib_results/ifeature15.tsv')
    df16 = pd.read_table ('biolib_results/ifeature16.tsv')
    df17 = pd.read_table ('biolib_results/ifeature17.tsv')
    df18 = pd.read_table ('biolib_results/ifeature18.tsv')
    df19 = pd.read_table ('biolib_results/ifeature19.tsv')
    df20 = pd.read_table ('biolib_results/ifeature20.tsv')
    df21 = pd.read_table ('biolib_results/ifeature21.tsv')

    MERGED_DATAFRAMES = []
    df_all = [df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16, df17, df18, df19, df20, df21]
    con1 = pd.concat(df_all, axis="columns")
    con1 = con1.loc[:, ~con1.columns.duplicated()] 
    MERGED_DATAFRAMES.append(con1)

    completed_steps += 1
    progress=int((completed_steps / total_steps) * 100)
    progress_callback("Step 3: Creating Dataset CSV", 4)

    con1.to_csv ('eskape_ifeature.csv', index=False)

    df_final = pd.read_csv('eskape_ifeature.csv')
    columns_to_keep = ['#','g3.g4.g3', 'AN.gap0', 'EM.gap0', 'LV.gap0', 'SV.gap0', 'QG.gap1', 'TI.gap1', 'FM.gap2', 'GQ.gap2', 'WC.gap2', 'RY.gap3', 'GD.gap4', 'YG.gap4', 'GS.gap5', 'PF.gap5', 'TR.gap5', 'TV.gap5', 'VI.gap5', 'ALF', 'AMP', 'APR', 'AVT', 'EFR', 'ELA', 'EYA', 'FKR', 'GGQ', 'GLR', 'GVP', 'HKG', 'IDD', 'IDR', 'IYQ', 'KLQ', 'KLS', 'LEK', 'NLT', 'NRY', 'NSG', 'PPV', 'RSL', 'RSS', 'SGE', 'SGH', 'TSV', 'VAG']
    final_df = df_final[columns_to_keep]
    final_df = final_df.rename(columns={'#': 'Protein_ID'})
    final_df.to_csv('ifeatures.csv', index=False)

    #features_file for users
    script_dir = sys.path[0]
    results_folder_path = os.path.join(script_dir, "Analysis_Results")

    #physicochemical features
    physico_filename = 'physiochemical_parameters.csv'
    physico_result_path = os.path.join(results_folder_path, physico_filename)
    con1.to_csv(physico_result_path, index=False)   

    fasta_file_path = "eskapeml/sequences.fasta"

    sequences = SeqIO.to_dict(SeqIO.parse(fasta_file_path, "fasta"))

    for header, sequence_record in sequences.items():
        sequence = str(sequence_record.seq)
        header_list.append(header)
 
        progress_callback("Sending POST request to MHC class I peptide binding prediction tool", 5)
        # Send POST request to MHC class I peptide binding prediction tool:
        mhci_res1 = iedb.query_mhci_binding(method="recommended", sequence= sequence, allele="HLA-A*02:01", length="9")
        mhci_res2 = iedb.query_mhci_binding(method="recommended", sequence= sequence, allele="HLA-A*01:02", length="9")

        mhci_res = pd.concat([mhci_res1, mhci_res2], axis=0)
        time.sleep(3)

    for header, sequence_record in sequences.items():
        sequence = str(sequence_record.seq)
        header_list.append(header)
    
        progress_callback("Sending POST request to MHC class II peptide binding prediction tool", 6)
        # Send POST request to MHC class II peptide binding prediction tool:
        mhcii_res1 = iedb.query_mhcii_binding(method="recommended", sequence= sequence, allele="HLA-DRB1*01:01", length=None)
        mhcii_res2 = iedb.query_mhcii_binding(method="recommended", sequence= sequence, allele="HLA-DRB1*01:02", length=None)
        
        mhcii_res = pd.concat([mhcii_res1, mhcii_res2], axis=0)
        time.sleep(3)

    for header, sequence_record in sequences.items():
        sequence = str(sequence_record.seq)
        header_list.append(header)
           
        progress_callback("Sending POST request to B-cell epitope prediction tool", 7)
        # Send POST request to B-cell epitope prediction tool:
        bcell_res = iedb.query_bcell_epitope(method="Bepipred", sequence= sequence, window_size=9)
        completed_steps += 1
        progress=int((completed_steps / total_steps) * 100)
        progress_callback("Step 4: Epitope Analysis Completed", 8)
        time.sleep(3)

    for header, sequence_record in sequences.items():
        sequence = str(sequence_record.seq)
        header_list.append(header)

        progress_callback("Sending POST request to adhesion probability prediction tool", 9)
        # Send POST request to adhesion probability prediction tool:
        sprob_res = iedb.query_bcell_epitope(method="Emini", sequence= sequence, window_size=9)
        completed_steps += 1
        progress=int((completed_steps / total_steps) * 100)
        progress_callback("Step 4: Adhesion Probability Analysis Completed", 10)
        time.sleep(3)

    for header, sequence_record in sequences.items():
        sequence = str(sequence_record.seq)
        header_list.append(header)

        progress_callback("Sending POST request to antigenicity prediction tool", 11)
        # Send POST request to antigenicity prediction tool:
        antigenicity_1 = iedb.query_bcell_epitope(method="Kolaskar-Tongaonkar", sequence= sequence, window_size=9)
        completed_steps += 1
        progress=int((completed_steps / total_steps) * 100)
        progress_callback("Step 5: Antigenicity Prediction Completed", 12)
        time.sleep(3)

        mhci_res['score'] = pd.to_numeric(mhci_res['score'], errors='coerce')
        mhci_res['percentile_rank'] = pd.to_numeric(mhci_res['percentile_rank'], errors='coerce')
            
        mhci_res_filtered = mhci_res[mhci_res["percentile_rank"] <= 0.5]
            
        if not mhci_res_filtered.empty:
            df_score_mhci = mhci_res_filtered["score"].mean()
            df_rank_mhci = mhci_res_filtered["percentile_rank"].mean()
            mhci_epitopes = mhci_res_filtered['peptide'].values
        else:
            df_score_mhci= 0
            df_rank_mhci= 0 

        mhcii_res['score'] = pd.to_numeric(mhcii_res['score'], errors='coerce')
        mhcii_res['rank'] = pd.to_numeric(mhcii_res['rank'], errors='coerce')
        
        mhcii_res_filtered = mhcii_res[mhcii_res["rank"] <= 1]

        if not mhcii_res_filtered.empty:
            rank_mhcii = mhcii_res_filtered["rank"].mean()
            score_mhcii = mhcii_res_filtered["score"].mean()
            mhcii_epitopes = mhcii_res_filtered['peptide'].values      
        else:
            rank_mhcii= 0
            score_mhcii= 0 

        bcell_res['Score'] = pd.to_numeric(bcell_res['Score'], errors='coerce')
        df_bcells_final = bcell_res["Score"].mean()

        sprob_res['Score'] = pd.to_numeric(sprob_res['Score'], errors='coerce')
        df_sprob_final = sprob_res["Score"].mean()

        antigenicity_1['Score'] = pd.to_numeric(antigenicity_1['Score'], errors='coerce')
        antigenicity_final = antigenicity_1["Score"].mean()

        progress_callback("Adding Scores to Dataset", 13)
        add_scores2 = pd.read_csv('ifeatures.csv')
        mhci_prob_scores.append(df_score_mhci)
        mhci_ranks.append(df_rank_mhci)
        mhcii_ranks.append(rank_mhcii)
        mhcii_scores.append(score_mhcii)
        b_cells_prob_scores.append(df_bcells_final)
        surface_probs.append(df_sprob_final)
        antigenicities.append(antigenicity_final)

        result_df = pd.DataFrame({
            #'Protein_ID': protein_ids,
            'antigenicity_1': antigenicities,
            'b_cells_probability_score': b_cells_prob_scores,
            'mhci_probability_score': mhci_prob_scores,
            'mhci_rank': mhci_ranks,
            'mhcii_rank': mhcii_ranks,
            'mhcii_score': mhcii_scores,
            'surface_probability': surface_probs,
        })
        result_df = result_df.round(6)
        merged_df = pd.concat([add_scores2, result_df], axis=1)
        merged_df.to_csv('ifeatures_updated.csv', index=False)

        completed_steps += 1
        progress=int((completed_steps / total_steps) * 100)
        progress_callback("Step 6: Scores added to Dataset", 14)

        
    progress_callback("Sending POST request to Signal peptides prediction tool", 15)
    # Send POST request to signal peptides prediction tool:
    signalp_6 = biolib.load('DTU/SignalP_6')
    signalp6_job = signalp_6.cli(args='--fastafile eskapeml/sequences.fasta --output_dir output')
    signalp6_job.save_files('biolib_results/signalP')

    script_dir = sys.path[0]
    signalP_folder_path = os.path.join(script_dir, "biolib_results/signalP")

    completed_steps += 1
    progress=int((completed_steps / total_steps) * 100)
    progress_callback("Step 7: Signal Peptide Analysis Completed", 16)

    sp_table_path = os.path.join(signalP_folder_path, "prediction_results.txt")
    df = pd.read_table(sp_table_path, sep="\t", header=None, skiprows=1)
    df.columns = ["ID", "Prediction", "OTHER", "SP(Sec/SPI)", "LIPO(Sec/SPII)", "TAT(Tat/SPI)", "TATLIPO(Sec/SPII)", "PILIN(Sec/SPIII)", "CS Position"]
    df.drop(index=df.index[0], axis=1, inplace=True)

    SP = df['SP(Sec/SPI)'].values
    LIPO = df['LIPO(Sec/SPII)'].values
    TAT = df['TAT(Tat/SPI)'].values
    TATLIPO = df['TATLIPO(Sec/SPII)'].values
    PILIN = df['PILIN(Sec/SPIII)'].values
    OTHER = df['OTHER'].values
    
    #features_file for users
    script_dir = sys.path[0]
    results_folder_path = os.path.join(script_dir, "Analysis_Results")

    threshold_signalp = [
                'SignalP result threshold: If a value of any column is 1 or close to 1, the protein belongs to that category',
            ]
    threshold_signalp_df = pd.DataFrame({'Thresholds': threshold_signalp})

    signalp_data = {
    'signal_peptide_SP': SP,
    'signal_peptide_LIPO': LIPO,
    'signal_peptide_TATLIPO': TATLIPO,
    'signal_peptide_TAT': TAT,
    'signal_peptide_PILIN': PILIN,
    'signal_peptide_OTHERS': OTHER  
    }

    df_signalp = pd.DataFrame(signalp_data)

    biological_filename = 'biological_properties.csv'

    biological_df =  pd.DataFrame({
        'Protein_ID': header,
        'antigenicity_1': antigenicities,
        'b_cells_probability_score': b_cells_prob_scores,
        'mhci_probability_score': mhci_prob_scores,
        'mhci_rank': mhci_ranks,
        'mhcii_rank': mhcii_ranks,
        'mhcii_score': mhcii_scores,
        'surface_probability': surface_probs,
    })

    biological_df = biological_df.round(6)


    threshold_adhesion = [
        'Threshold for adhesion probability = 1 Or greater than 1',
    ]
    threshold_adesion_df = pd.DataFrame({'Thresholds': threshold_adhesion})

    
    df_mhci_epitopes = pd.DataFrame({
                                        'Epitope_Type': ['mhci'] * len(mhci_epitopes),
                                        'Epitope_Sequence': mhci_epitopes})
    
    df_mhcii_epitopes = pd.DataFrame({
                                        'Epitope_Type': ['mhcii'] * len(mhcii_epitopes),
                                        'Epitope_Sequence': mhcii_epitopes})
    
    comments = [
        'Threshold for strong "MHC Class I" binder: % Rank = Lower than 0.5',
        'Threshold for strong "MHC Class II" binder: % Rank = Lower than 1',
        'Threshold for strong "B cell epitope": Score = Greater than 0.5',
        'Antigenicity = Value close to 1 or Greater than 1 indicates high antigenicity'
    ]
    
    biological_result_path = os.path.join(results_folder_path, biological_filename)
    
    comments_df = pd.DataFrame({'Thresholds': comments})

    df_biological = pd.concat([biological_df, df_signalp, df_mhci_epitopes, df_mhcii_epitopes, comments_df, threshold_signalp_df], ignore_index=True)
        
    df_biological.to_csv(biological_result_path, index=False)

    progress_callback("Step 8: Adding SignalP Scores to Dataset", 17)

    add_scores = pd.read_csv('ifeatures_updated.csv')
    add_scores["signal_peptide_SP"] = SP
    add_scores["signal_peptide_LIPO"] = LIPO
    add_scores["signal_peptide_TAT"] = TAT
    add_scores["signal_peptide_TATLIPO"] = TATLIPO
    add_scores["signal_peptide_PILIN"] = PILIN
    add_scores["signal_peptide_OTHER"] = OTHER

    add_scores.to_csv('finalfeatures_updated.csv', index=False)


    progress_callback("Step 9: DataSet Created Successfully", 18)

def upload_sequence(request):
    if request.method == "POST":
        sequence = request.POST.get("sequence")
        #sequence = sequence.replace("\n", "")
        file = request.FILES.get("file")

        try:
            if sequence:
                script_dir = os.path.dirname(os.path.abspath(__file__))
                file_path = os.path.join(script_dir, "sequences.fasta")

                with open(file_path, "w") as f:
                    f.write(sequence)
            elif file:
                script_dir = os.path.dirname(os.path.abspath(__file__))
                file_path = os.path.join(script_dir, "sequences.fasta")
                with open(file_path, "wb") as f:
                    for chunk in file.chunks():
                        f.write(chunk)

            calculate_features(request, file_path, progress_callback)
            

        except Exception as e:
            print(f"Error analyzing sequence: {e}")
            return JsonResponse({"status": "Error analyzing sequence"}, status=500)

        return redirect("results")
    
    # Logger and Progress Bar
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    script_dir = os.path.dirname(os.path.abspath(__file__))
    log_file_path = os.path.join(script_dir, "vacsol_ml.log")
    file_handler = logging.FileHandler(log_file_path)
    file_handler.setLevel(logging.INFO)
    formatter = ProgressBarFormatter('%(message)s')
    file_handler.setFormatter(formatter)
    logger.addHandler(file_handler)

    with open('vacsol_ml.log', 'r') as log_file:
        logs = log_file.readlines()

    with open('vacsol_ml.log', 'w') as log_file:
        log_file.truncate(0)
        
    with open('eskapeml/vacsol_ml.log', 'w') as log_file2:
        log_file2.truncate(0)

    return render(request, 'upload.html', {'logs': logs})

class ProgressBarFormatter(logging.Formatter):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.completed_steps = 0
        self.total_steps = 10
        self.progress_bar_added = True

    def format(self, record):
        if hasattr(record, 'progress'):
            self.completed_steps = record.progress

        progress = int((self.completed_steps / self.total_steps) * 100)

        progress_bar = f"[{'#' * (progress // 10)}{'-' * ((100 - progress) // 10)}] {progress}%"

        if not self.progress_bar_added:
            record.msg = f"{progress_bar}\n{record.msg}"
            self.progress_bar_added = True
        else:
            record.msg = f"\n{record.msg}"

        return super().format(record)
        
# Results
def get_results(request):
    try:
        model_path = 'eskapeml/ML_model.pkl'

        with open(model_path, "rb") as file:
            model = joblib.load(file)

        final = pd.read_csv('finalfeatures_updated.csv')

        Protein_ID = final['Protein_ID']
        finalDF = final.drop(labels=["Protein_ID"], axis=1)

        output_pred = model.predict(finalDF)

        output_proba = model.predict_proba(finalDF)

        proba_class_1 = output_proba[:, 1]
        proba_class_0 = output_proba[:, 0]


        results = pd.DataFrame({"Protein_ID": final['Protein_ID'], "Prediction": output_pred, "Probability_Class_1": proba_class_1, "Probability_Class_0": proba_class_0,
                                "antigenicity_1": final["antigenicity_1"], "b_cells_probability_score": final["b_cells_probability_score"],
                                "mhci_probability_score": final["mhci_probability_score"], "mhci_rank": final["mhci_rank"], "mhcii_rank": final["mhcii_rank"], "mhcii_score": final["mhcii_score"],
                                 "surface_probability": final["surface_probability"], "signal_peptide_SP": final["signal_peptide_SP"], "signal_peptide_LIPO": final["signal_peptide_LIPO"], "signal_peptide_TAT": final["signal_peptide_TAT"], "signal_peptide_TATLIPO": final["signal_peptide_TATLIPO"], "signal_peptide_PILIN": final["signal_peptide_PILIN"],
                                })

        results_data = [
            {"Protein_ID": row["Protein_ID"], "Prediction": int(row["Prediction"]), "Probability_Class_1": row["Probability_Class_1"], "Probability_Class_0": row["Probability_Class_0"],
             "antigenicity_1": row["antigenicity_1"], "b_cells_probability_score": row["b_cells_probability_score"],
             "mhci_probability_score": row["mhci_probability_score"], "mhci_rank": row["mhci_rank"], "mhcii_rank": row["mhcii_rank"], "mhcii_score": row["mhcii_score"],
             "surface_probability": row["surface_probability"], "signal_peptide_SP": row["signal_peptide_SP"], "signal_peptide_LIPO": row["signal_peptide_LIPO"], "signal_peptide_TAT": row["signal_peptide_TAT"], "signal_peptide_TATLIPO": row["signal_peptide_TATLIPO"], "signal_peptide_PILIN": row["signal_peptide_PILIN"],
             }
            for _, row in results.iterrows()
        ]

        return render(request, "results.html", {"results": results_data})

    except Exception as e:
        print(f"Error analyzing sequence: {e}")
        return JsonResponse({"status": "Error analyzing sequence"}, status=500)
