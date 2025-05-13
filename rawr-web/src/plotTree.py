#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# This file is used for MSA and tree support estimation.

import ete3
import os
import argparse


def plotTreeWithSupport(treeFile, treeFigFile):
    ts = ete3.TreeStyle()
    ts.show_leaf_name = True
    ts.show_branch_length = True
    ts.show_branch_support = True
    if os.path.exists(treeFile):
        with open(treeFile, "r") as inf:
            try:
                tree = ete3.Tree(inf.readline())
                tree.render(treeFigFile, w=183, units="mm", tree_style=ts)
                print("finished.")
            except BaseException:
                print("No valid tree found in " + treeFile + ".")

def main():
    parser = argparse.ArgumentParser(
        description='Parameters for calculating support for MSA.')
    parser.add_argument(
        '--treeFile',
        '-t',
        help='input tree file in .newick format.',
        required=True)
    parser.add_argument(
        '--figFile',
        '-f',
        help='output figure file.',
        required=True)
    args = parser.parse_args()

    treeFile = args.treeFile
    figFile = args.figFile
    plotTreeWithSupport(treeFile, figFile)

    return


if __name__ == '__main__':
    main()
