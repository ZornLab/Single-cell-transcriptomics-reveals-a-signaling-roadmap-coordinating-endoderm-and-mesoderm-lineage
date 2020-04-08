import sys
import os
path = os.getcwd()
os.chdir(path)
os.system('python2.7 -m SimpleHTTPServer 8000 &')
#####ps -fA | grep python