import os
import exceptions
import re
from collections import OrderedDict

import numpy as np


class Namelist(object):

    def __init__(self, filename=None):
        self.top     = OrderedDict()
        self.valdict = OrderedDict()
        self.groups  = []

        if filename:
            self.read(filename)
            
    def read(self, filename, append=False):
        groupname = "top"
        for l in open(filename):
            cmd = l.split('!')[0].strip()
            if cmd:
                if ("&" in cmd) & (len(cmd.strip().split("&")[0]) == 0):
                    groupname =  cmd.strip().strip('&')
                    self.addgroup(groupname)
                elif ("$" in cmd) & (len(cmd.strip().split("$")[0]) == 0):
                    groupname =  cmd.strip().strip('$')
                    self.addgroup(groupname)

                elif "=" in cmd:
                    var,arg = cmd.split("=", 1)
                    self.addval(var, arg, groupname)
    
    def duck (self, instr):
        if "," in instr:
            try:
                return [float(s) for s in instr.split(',')]
            except:
                return instr.strip("'")
        try:
            return int(instr)
        except exceptions.ValueError:
            try:
                return float(instr)
            except exceptions.ValueError:
                return instr.strip("'")

    def addgroup(self, groupname):
        groupname = groupname.upper()
        if not hasattr(self, groupname):
            setattr(self, groupname, OrderedDict())
            self.groups.append(groupname.upper())
            
    def addval(self, var, arg, groupname="top"):
        """Add value to namelist"""
        arg = arg.strip().rstrip(',').rstrip('/')
        var = var.strip().lower()
        if "(" and ")" in var:
            ind = int(var[var.index("(") + 1:var.rindex(")")])-1
            listvar = var.split('(')[0].strip()
            if not listvar in self.__dict__:
                setattr(self,listvar,defaultlist())
            self.__dict__[listvar][ind] = self.duck(arg)
            self.__dict__[groupname][listvar] = self.__dict__[listvar]
            self.valdict[listvar] = self.__dict__[listvar]

        else:
            self.__dict__[groupname][var] = self.duck(arg)
            setattr(self, var, self.duck(arg))
            self.valdict[var] = self.duck(arg)
                    
    def iternml(self):
        "Create a iterobject to generate a nml file"""
        for group in self.groups:
            yield "&%s" % group
            for key,val in self.__dict__[group].items():
                if type(val) is defaultlist:
                    if len(val) >1:
                        val = ",".join(str(e) for e in val)
                    else:
                        val = val[0]
                elif type(val) is str:
                    val = "'%s'" % val
    
                yield '{0: >16} = {1: >6},'.format(key,val)
            yield "/"

    def pprint(self):
        for l in self.iternml():
            print l

    def write(self,filename):

        with open(filename, 'w') as fH:
            for l in self.iternml():
                fH.write(l + "\n")


    def __repr__(self):
        self.pprint()
        return ""
            
class defaultlist(list):

    def __setitem__(self, index, value):
      size = len(self)
      if index >= size:
         self.extend(0 for _ in range(size, index + 1))

      list.__setitem__(self, index, value)



def upgrade(templatename, infilename, outfilename, explist=[]):

    nlt = Namelist(templatename)
    if type(infilename) is list:
        nli = Namelist(infilename[0])
        for name in infilename[1:]:
            nli.read(name)
    else:
        nli = Namelist(infilename)

        
    for group in nlt.groups:
        for key in nlt.__dict__[group].keys():
            if hasattr(nli, key) & (not key in explist):
                arg = getattr(nli, key)
                nlt.__dict__[group][key] = arg
                setattr(nlt, key, arg)
                nlt.valdict[key] = arg
    nlt.write(outfilename)


def tracmass_upgrade_all_projects():

    for d in os.walk('./'):
        stamp = "%s/%s%s.in"
        if os.path.isfile(stamp % (d[0], d[0], "_run")):
            print d[0]
            infilenames = [stamp % (d[0], d[0], "_run"),
                           stamp % (d[0], d[0], "_grid")]
            outfilename = stamp % (d[0], d[0], "")
            upgrade("template.in", infilenames, outfilename, ['gridvernum',])
