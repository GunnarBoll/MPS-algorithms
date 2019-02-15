import pathlib
import os

def cwd_store(dname, fname, data):
    
    run_nr = 1
    path_direc = (os.getcwd() + "/" + dname)
    while True:
        if run_nr > 10:
            break
        try:
            direc = path_direc + "_" + str(run_nr) + "/"
            pathlib.Path(direc).mkdir(parents=True, exist_ok=True)
            with open(direc+fname, "x") as fw:
                for mat in data:
                    fw.write(str(mat) + "\n")
            break
        except FileExistsError:
            run_nr += 1
    
    return