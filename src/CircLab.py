#!/usr/bin/env python
# -*- coding: UTF-8 -*-
'''
@Project ：pcircom
@File    ：CircLab.py
@Author  ：xm
@Date    ：2023/12/4 下午5:15
'''
#!/usr/bin/env python
# -*- coding: utf-8 -*-

'''
Usage: ToolName (<command> | --help | --version)
'''

from docopt import docopt
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.realpath(__file__)))

#from version import __version__

__author__ = 'xm (your.email@example.com)'

help_doc = r'''
		  _                     __            __        
		 (_)                   [  |          [  |       
	 .---.   __    _ .--.   .---.   | |   ,--.    | |.--.   
	/ /'`\] [  |  [ `/'`\] / /'`\]  | |  `'_\ :   | '/'`\ \ 
	| \__.   | |   | |     | \__.   | |  // | |,  |  \__/ | 
	'.___.' [___] [___]    '.___.' [___] \'-;__/ [__;.__.'  

   circlab：A comprehensive software tool for circRNA detection, 
            function prediction, functional enrichment, and visualization.
                                                              
   Usage: GreenCircLab <command> [options] - Choose a module to execute

   Module:
       detection       - Detect circRNAs from sequencing data using various prediction tools and generate unified results.
       function        - Predict potential interactions between circRNAs and miRNAs/mRNAs and the IRES element in circRNA.
       annotation      - Annotate circRNAs by parental genes and predicted mRNA targets to functional enrichment (GO/KEGG).
       visualization   - Generate visual plots for circRNA classification, expression distribution, and genomic context.
'''


def main():
    # parse command
    command_log = 'GreenCircLab parameters: ' + ' '.join(sys.argv)
    if len(sys.argv) == 1:
        sys.exit(help_doc)
    # elif sys.argv[1] == '--version' or sys.argv[1] == '-v':
    #     sys.exit(__version__)
    elif sys.argv[1] == 'detection':
        print("start detection")
        import detection
        detection.main(sys.argv[2:])
    elif sys.argv[1] == 'function':
        import function
        function.main(sys.argv[2:])
    elif sys.argv[1] == 'annotation':
        import annotation
        annotation.main(sys.argv[2:])
    elif sys.argv[1] == 'visual':
        import visual
        visual.main(sys.argv[2:])
    else:
        sys.exit(help_doc)


if __name__ == '__main__':
    main()
