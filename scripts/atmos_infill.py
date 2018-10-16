# -*- coding: utf-8 -*-
'''
Infill CMIP5 atmospheric fields over land

Created on Mon Dec  7 08:39:36 2015

@author: James Harle
'''

import sys, getopt
import infill

def main():
    """ Main function which checks the command line parameters
        before calling the main function"""

    atm_file = ''
    msk_file = ''
    var_name = ''
    int_meth = ''
    try:
        opts, dummy_args = getopt.getopt(sys.argv[1:], "hf:v:m:i:", 
                                         ["help","file=","variable","mask=","method="])
    except getopt.GetoptError:
        print "usage: atmos_infill -f <input.nc> -v <variable> -m <mask_file.nc -i <interpolation method>"
        sys.exit(2)

    for opt, arg in opts:
        if opt == "-h":
            print "usage: atmos_infill -f <input.nc> -v <variable> -m <mask_file.nc -i <interpolation method>"
            print "       -f <input.nc>       input file to use"
            print "       -v <variable>       variable to extract"
            print "       -m <mask_file.nc>   mask file to use"
            print "       -i <interpolation method>   [Rbf|flood]"
            sys.exit()
        elif opt in ("-f", "--file"):
            atm_file = arg
        elif opt in ("-m", "--mask"):
            msk_file = arg
        elif opt in ("-v", "--variable"):
            var_name = arg
        elif opt in ("-i", "--method"):
            int_meth = arg

    if atm_file == "":
        print "usage: atmos_infill -f <input.nc> -v <variable> -m <mask_file.nc -i <interpolation method>"
        sys.exit(2)
    if msk_file == "":
        print "usage: atmos_infill -f <input.nc> -v <variable> -m <mask_file.nc -i <interporlation method>"
        sys.exit(2)
    if var_name == "":
        print "usage: atmos_infill -f <input.nc> -v <variable> -m <mask_file.nc -i <interpolation method>"
        sys.exit(2)
    if int_meth == "":
        print "usage: atmos_infill -f <input.nc> -v <variable> -m <mask_file.nc -i <interpolation method>"
        sys.exit(2)
    
    if int_meth.lower() == "rbf":
        infill.infill_Rbf(atm_file,var_name,msk_file)
    elif int_meth.lower() == "flood":
        infill.infill_flood(atm_file,var_name,msk_file)
    else:
        print "interpolation method is either Rbf or flood"
        
if __name__ == "__main__":
    main()
