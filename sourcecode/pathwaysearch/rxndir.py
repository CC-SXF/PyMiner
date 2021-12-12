# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 11:59:30 2020

@author: CC-SXF
"""

import re
import json

class RxnDir():
    """
    """
    def __init__(self):
        """ """
        pass


    @classmethod
    def _getRxnDirection(cls, rxn_equ):
        """ """
        # return "=>"
        # return "<="
        # return "<=>"
        pattern_equ = re.compile("( <=> | => | <= | = )")
        rxn_dir = re.search(pattern_equ, rxn_equ).group(1)
        rxn_dir = rxn_dir.strip()
        if rxn_dir == "=":
            rxn_dir = "<=>"
        rxn_dir = json.dumps(rxn_dir)
        return rxn_dir


    @classmethod
    def demoFunc(cls,):
        """ """
        pass


if __name__ == "__main__":
    """ """
    pass





