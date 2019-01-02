import os, sys, numpy, tabulate
import pandas as pd
from sklearn.metrics import precision_recall_curve, roc_curve
from sklearn import metrics
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()
from snakemake.utils import R

SECRETS = '/home/evansj/me/.secrets/'
sys.path.append(SECRETS)
from pass_wd import *

DONE = TWILIO_PRE + "--data-urlencode 'Body=DONE' " + TWILIO_POST
FAIL = TWILIO_PRE + "--data-urlencode 'Body=FAIL' " + TWILIO_POST

p = os.getcwd()
if 'src' in p:
    PWD = p.split('src/')[0]
else:
    PWD = p + '/'

FILES = PWD + 'docs/'
SCRIPTS = PWD + 'src/scripts/'
DATA = PWD + 'data/'
LOG = PWD + 'log/'
CONFIG = PWD + 'configs/'
PLOTS = PWD + 'writeup/plots/'
ENVS = PWD + 'envs/'

# containers
PLINK = "docker://quay.research.chop.edu/evansj/plink-docker:181012"
PLINK2 = "docker://quay.research.chop.edu/evansj/plink2-docker:11Feb"
PRSICE = "docker://quay.research.chop.edu/evansj/prsice-docker:2.1.4"

G = ('early', 'all', 'late', 'ibd_all')
