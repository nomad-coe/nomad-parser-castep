# CASTEP Parser

This is the parser for [CASTEP](http://www.castep.org/).
It is part of the [NOMAD Laboratory](http://nomad-lab.eu).
The official version lives at

    git@gitlab.mpcdf.mpg.de:nomad-lab/parser-castep.git

you can browse it at

    https://gitlab.mpcdf.mpg.de/nomad-lab/parser-castep

It relies on having the nomad-meta-info and the python common repositories one level higher.
The simplest way to have this is to check out nomad-lab-base recursively:

    git clone --recursive git@gitlab.mpcdf.mpg.de:nomad-lab/nomad-lab-base.git

then this will be in parsers/castep.


# TEST OUTPUT FILES 

Few output files where to test the parser are provided in the directory test examples.

        FILE NAME     |              FILE DESCRIPTION
    __________________|___________________________________________________
    "Si2.castep_v_1" --> Single Point Calculation (minimum verbosity)
    "Si2.castep_v_2" --> Single Point Calculation (medium verbosity)
    "Si2.castep_v_3" --> Single Point Calculation (maximum verbosity)
    
    "Si2.castep_b_v_1" --> Band Structure Calculation (minimum verbosity)
    "Si2.castep_b_v_2" --> Band Structure Calculation (medium verbosity)
    "Si2.castep_b_v_3" --> Band Structure Calculation (maximum verbosity)