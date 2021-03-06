# PyMiner
PyMiner is a practical and effective tool for metabolic pathway design that can meet a variety of pathway design requirements, including the pathways from initial substrates to target product, and the exogenous pathways of one specific chassis microorganism just given the target product. In addition, PyMiner can search for metabolic pathways within a given length, metabolic pathways with a specific length , and the shortest metabolic pathways.

<img src="https://github.com/CC-SXF/PyMiner/blob/main/PyMiner.jpg" width="800" />

# Anaconda environment

+ Windows:

  https://repo.anaconda.com/archive/Anaconda3-2019.07-Windows-x86_64.exe

  -> Add Anaconda to the system PATH environment variable.
+ MacOS:

  https://repo.anaconda.com/archive/Anaconda3-2019.07-MacOSX-x86_64.pkg

PyQt5, RDKit (https://www.rdkit.org/) and COBRApy (http://opencobra.sourceforge.net/) were used to build the metabolic pathway design tool PyMiner.

+ PyQt5 (satisfied after installing Anaconda3)
  
  https://pypi.org/project/PyQt5/5.9.2/
```
  # # pip install PyQt5==5.9.2
```

+ NumPy (so as to be compatible with COBRApy)
  
  https://numpy.org/install/
```
  pip install numpy==1.17.3
```

+ RDKit
  
  https://anaconda.org/rdkit/rdkit
```
  conda install -c rdkit rdkit=2020.03.2
```

+ COBRApy
  
  https://pypi.org/project/cobra/0.22.0/
```
  pip install cobra==0.22.0
```

# Usage
```
PyMiner.py
```

# References
Song XF, Dong MY, Liu M. PyMiner: A Method for Metabolic Pathway Design Based on the Uniform Similarity of Substrate-Product Pairs and Conditional Search. PLoS ONE. 2022; 17(4):e0266783. https://doi.org/10.1371/journal.pone.0266783
