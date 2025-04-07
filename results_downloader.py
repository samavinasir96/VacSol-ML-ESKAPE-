import os
import sys
import zipfile
from django.http import HttpResponse
from django.shortcuts import render


def download_csv(request):
    if request.method == 'POST':
        selected_files = request.POST.getlist('csv_files')
        if selected_files:
            # Set the response content type as a zip file
            response = HttpResponse(content_type='application/zip')
            response['Content-Disposition'] = 'attachment; filename="analysis_results.zip"'
            # Create a zip file to store the selected CSV files
            with zipfile.ZipFile(response, 'w') as zip_file:
                # Path to the directory containing the CSV files
                script_dir = sys.path[0]
                directory = os.path.join(script_dir, 'Analysis_Results')
                for file_name in selected_files:
                    # Path to the CSV file
                    file_path = os.path.join(directory, file_name)
                    # Add the CSV file to the zip file
                    zip_file.write(file_path, os.path.basename(file_path))
            return response
    else:
        return render(request, 'results.html')
