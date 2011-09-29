import exceptions

import numpy as np
import pylab as pl

def duck (s):
    try:
        return int(s)
    except exceptions.ValueError:
        try:
            return float(s)
        except exceptions.ValueError:
            return s.strip("'")

def parse(filename):
    class struct:
        pass
    struct.top = {}
    curstr = 'top'
    for l in open(filename):
        cmd = l.split('!')[0].strip()
        if cmd:
            if "&" in cmd:
                dlm =  cmd.strip().strip('&')
                struct.__dict__[dlm] = {}
                curstr = dlm
            elif "=" in cmd:
                var,arg = cmd.split("=")
                arg = arg.strip().rstrip(',').rstrip('/')
                var = var.strip()
                struct.__dict__[curstr][var] = duck(arg)
                print var,type(arg)
    return struct
