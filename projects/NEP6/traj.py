import pytraj

import pdb; pdb.set_trace()
tr = pytraj.Trm('NEP6', 'NEP6', 'bering')
#tr.load()                      #Load the last generated file
#tr.load(jd=731583)             #Load the file starting at 1234
tr.load(filename="NEP6_test00731582_run.asc") #Load the file named file.bin

tr.scatter(jd=tr.jdvec[0])
