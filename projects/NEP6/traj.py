import pytraj

tr = pytraj.Trm('NEP6', 'NEP6')
tr.load()                      #Load the last generated file
tr.load(jd=731583)             #Load the file starting at 1234
tr.load (filename="NEP6_test00731582_run.asc") #Load the file named file.bin   
