import rpy2
from rpy2.robjects.packages import importr

packnames = ('dendextend', 'pvclust', 'grDevices', 'ape')

utils = importr('utils')
utils.chooseCRANmirror(ind=1)
for package in packnames:
    if package not in rpy2.robjects.r['installed.packages']():
        utils.install_packages(package)