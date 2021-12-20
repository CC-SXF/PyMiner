# PyMiner
PyMiner is a practical and effective tool for metabolic pathway design that can meet a variety of pathway design requirements, including the pathways from initial substrates to target product, and the exogenous pathways of one specific chassis microorganism just given the target product. In addition, PyMiner can search for metabolic pathways within a given length, metabolic pathways with a specific length , and the shortest metabolic pathways.

# Anaconda environment

+ Windows:
  https://repo.anaconda.com/archive/Anaconda3-2019.07-Windows-x86_64.exe
  -> Add Anaconda to the system PATH environment variable.
+ MacOS:
  https://repo.anaconda.com/archive/Anaconda3-2019.07-MacOSX-x86_64.pkg

PyQt5, RDKit (https://www.rdkit.org/) and COBRApy (http://opencobra.sourceforge.net/) were used to build the metabolic pathway design tool PyMiner.

+ 1) PyQt5 (satisfied after installing Anaconda3)
  
  https://pypi.org/project/PyQt5/5.9.2/
```
# # pip install PyQt5==5.9.2
```

+ 2) numpy (so as to be compatible with COBRApy)
  
  https://numpy.org/install/
```
pip install numpy==1.17.3
```

+ 3) RDKit
  
  https://anaconda.org/rdkit/rdkit
```
conda install -c rdkit rdkit=2020.03.2
```

+ 4) COBRApy
  
  https://pypi.org/project/cobra/0.22.0/
```
pip install cobra==0.22.0
```

# Usage
```
PyMiner_1.0.py
```

# References
If you find this work or code useful, please cite this study. The citation will be updated soon.
