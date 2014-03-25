import readline, glob, os,sys
        

def complete(text, state):
    return (glob.glob(text+'*')+[None])[state]
readline.set_completer_delims(' \t\n;')
readline.parse_and_bind("tab: complete")
readline.set_completer(complete)

print("checking R in path..."), 
R = 0
for path in os.environ["PATH"].split(":"):
    if '/R/' in path : 
        print "[OK]"
        R += 1
if R == 0 :
    print "Could not find R in $PATH"
    sys.exit()

print("checking subprocess..."),
try:
   import subprocess
except ImportError:
   print >> sys.stderr,"Could not import the subprocess module"
   raise
print "[OK]"

print("checking RPy2..."),
try:
   from rpy2.robjects.packages import importr
except ImportError:
   print >> sys.stderr,"Could not import RPy2"
   raise
print "[OK]"


