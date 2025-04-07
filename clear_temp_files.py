# clear files after user downloaded them, to clear up the server

import os
import sys

from django.http import HttpResponseRedirect


def delete_files(request):
    script_dir = sys.path[0]
    folder_path = os.path.join(script_dir, 'biolib_results')
    folder_path2 = os.path.join(script_dir, 'Analysis_Results')
    folder_path3 = os.path.join(script_dir, 'biolib_results', 'signalP')

    # Delete all files in the folder
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        if os.path.isfile(file_path):
            os.remove(file_path)
    
    for filename in os.listdir(folder_path2):
        file_path2 = os.path.join(folder_path2, filename)
        if os.path.isfile(file_path2):
            os.remove(file_path2)
    
    for filename in os.listdir(folder_path3):
        file_path3 = os.path.join(folder_path3, filename)
        if os.path.isfile(file_path3):
            os.remove(file_path3)

    return HttpResponseRedirect('/')
