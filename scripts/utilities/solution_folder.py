# importing the latest folder from the solutions folder

import os
import sys

# The path to the solution folder
soln_path = os.path.join(os.path.dirname(os.path.split(sys.path[1])[0]), 'solutions')

# Get the full path to the solution folder_name
def getpath(folder_name) :

    return os.path.join(soln_path, folder_name)

# Returns the last modified time of the solution folder_name
def maxkeyfunc(folder_name) :

    return os.path.getmtime(getpath(folder_name))

# Get the latest solution created / modified
def getlatestfolder() :

    return getpath(max(os.listdir(soln_path), key=maxkeyfunc))