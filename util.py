import os

#generate all the directories needed for the given path
def generateDirectories(path):
	folders=path.split("/")
	curdir=""
	for folder in folders:
		curdir=os.path.join(curdir,folder)
		if not os.path.exists(curdir):
			os.mkdir(curdir)
def generateDirectoriesMult(paths):
	for path in paths:
		generateDirectories(path)