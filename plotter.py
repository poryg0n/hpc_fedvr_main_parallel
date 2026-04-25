import os 

with os.scandir('data/') as data_dir:
    for dir_ in data_dir:
        print(dir_)
