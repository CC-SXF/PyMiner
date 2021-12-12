# -*- coding: utf-8 -*-
"""
Created on Sat Dec 11 15:19:26 2021

@author: CC-SXF
"""

from os import path
import tarfile


def getMembers(members,  suffix='.rxn'):
    """
    """
    for tarinfo in members:
        if tarinfo.name.endswith(suffix):
            tarinfo.name = path.split(tarinfo.name)[-1]
            yield tarinfo


def unpackResource(compressed_file, target_path, suffix='.rxn'):
    """
    """
    tar = tarfile.open(compressed_file)
    members = getMembers(tar, suffix)
    tar.extractall(path=target_path, members=members)
    tar.close()


if __name__ == '__main__':
    """
    """
    unpackResource('./temp/KEGG_Debug.tar.gz', './temp/','.data')
    unpackResource('./temp/KEGG_Dmn.tar.gz', './temp/','.data')

    unpackResource('./temp/KndPad_Debug.tar.gz', './temp/','.data')
    unpackResource('./temp/KndPad_Dmn.tar.gz', './temp/','.data')

    unpackResource('./temp/MetaCyc_Debug.tar.gz', './temp/','.data')
    unpackResource('./temp/MetaCyc_Dmn.tar.gz', './temp/','.data')

    unpackResource('../KEGGAAMBCsRCsFiles_Ori_Half.tar.gz', '../KEGGAAMBCsRCsFiles_Ori_Half')

    unpackResource('../KndPadAAMBCsRCsFiles_Ori_Half_part_1.tar.gz', '../KndPadAAMBCsRCsFiles_Ori_Half')
    unpackResource('../KndPadAAMBCsRCsFiles_Ori_Half_part_2.tar.gz', '../KndPadAAMBCsRCsFiles_Ori_Half')
    unpackResource('../KndPadAAMBCsRCsFiles_Ori_Half_part_3.tar.gz', '../KndPadAAMBCsRCsFiles_Ori_Half')





