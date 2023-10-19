# :mortar_board: MSc_PGTKEP_IndividualRepo

This repository consists of the contents for my individual PGTKEP project. 

Below description illustrates in details. 

<!--
### :computer: Reference

The original source code is publicly available on:
https://github.com/expartools/ETSeq
-->

### :bookmark_tabs: Experiments

#### Replicate_ETSeq
> The original source code is publicly available on: https://github.com/expartools/ETSeq

:file_folder: Data files

* templateSeqLists.txt
  > Dataset of all template sequences from Qian et al.

* templateSeqLists_EXCEPT_77seq.txt
  > Dataset of template sequences without 77 templates
  
:computer: Replication
  
The `run.py` is the key file to run the computational tool called `ETSeq`. To run `run.py`, you need to prepare environment as below:
* Create a python virtual environment named .venv

  `python -m venv .venv`
  
* Activate your python virtual environment 

  Mac/Linux `source .venv/bin/activate`
    
  Windows `.venv/Scripts/activate.bat`
* Install `Orange` to run `run.py`
  
  `pip install PyQt5 PyQtWebEngine`

  `pip install orange3`

#### Reproduce_ETSeq
> The original source code is publicly available on: https://github.com/expartools/ETSeq

:file_folder: Data files

* templateSeqs_withonlyClass.xlsx
  > Combination of Class I and II template sequences (from Qian et al.)

#### Implement_LDpipeline
I've conducted the computational pipeline for  the development of a Point-of-Care (PoC) testing procedure for sexually transmitted infections (STIs) at Linear Diagnostics Limited. 
However, **due to security reason**, this repository is private.

Please contact me or Jean-Louis Duprey (Head of R&D at Linear Diagnostics) for the further questions. :speech_balloon:
