## To export MATLAB F0AM output structure into Python



Download the files in this GitHub repository as a directory:

**https://github.com/jerelbn/yamlmatlab**



This is based off of **\*this original package.https://code.google.com/archive/p/yamlmatlab/** If unavailable, there is also a GitHub repository of the same package at: **\*https://github.com/ewiger/yamlmatlab\***



Make sure to add the directory to your MATLAB path.



In order to use, insert into your MATLAB script:

WriteYaml(yaml\_save\_path,S);



where yaml\_save\_path is a string of your file's savepath \& name (so it should end in '.yaml') and S is the name of the MATLAB structure that you want to save as a yaml file.



## To read a YAML file in Python

Download this package (preferably in your conda/mamba environment):

**\*html link\***



To install it with micromamba, use:

**\*insert code here\***



To import the package, use

import yaml



even though the package is called pyyaml, it is referred to in the system as only "yaml" when importing.



To open the file, define the path for your yaml (I am using yaml\_file\_name) and use the following code:



with open(yaml\_file\_name, 'r') as f:

    yaml\_variable = yaml.full\_load(f)



where yaml\_variable is the variable I will use to hold my yaml data.



