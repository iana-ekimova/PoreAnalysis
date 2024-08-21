# Pore Analysis

## Installation

Clone repository and update submodules:
`git clone git@github.com:iana-ekimova/PoreAnalysis.git -b main`
`git submodule init`
`git submodule update`

Install dependencies:
`pip install requirements.txt`

Install pypore3d:
`cd 3rdparty/pypore3d`
`python setup.py build_ext`
`cd ../../`
`pip install 3rdparty/pypore3d`

### Troubleshooting on Mac:
There could be problems with a `libomp`, so it needed to run this command (after installation of libomp with brew):
`python setup.py build_ext --inplace --include-dirs=/opt/homebrew/Cellar/libomp/18.1.5/include`