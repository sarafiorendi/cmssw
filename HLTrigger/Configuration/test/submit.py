import os
import datetime
from optparse import OptionParser

parser = OptionParser()
parser.usage = '''
'''
parser.add_option("-Q"  , "--queue"   , dest = "queue"    ,  help = "choose queue. Available 1nh 8nh 1nd 2nd 1nw 2nw. Default is 1nd" , default = '1nd' )
parser.add_option("-N"  , "--njobs"   , dest = "njobs"    ,  help = "choose number of jobs.Default is -1"         , default = -1                             )
parser.add_option("-F"  , "--nfiles"  , dest = "nfiles"   ,  help = "choose number of files per job.Default is 1" , default =  1                        )
parser.add_option("-S"  , "--start"   , dest = "start"    ,  help = "choose starting file"                        , default =  0                             )
parser.add_option("-i"  , "--input"   , dest = "input"    ,  help = "input file list"                             , default =  "../fileList.txt"          )
parser.add_option("-c"  , "--cfg"     , dest = "cfg"      ,  help = "input cfg file"                              , default =  "../hlt_cfg.py"            )
parser.add_option("-d"  , "--outdir"  , dest = "outdir"   ,  help = "output dir on eos"                           , default =  "PFMuonHLT_forTrackComparison_veto0p01Charged_hitsCutAndPt0p5_noDZ")
parser.add_option("-t"  , "--hlt"     , dest = "hlt"      ,  help = "hlt path "                                   , default =  "HLT_PFIsoMu24_eta2p1"        )

(options,args) = parser.parse_args()	

njobs  = options.njobs
nfiles = options.nfiles  

njobs  = int(njobs)
nfiles = int(nfiles)




myTime = []
for i in list(datetime.datetime.now().timetuple())[:5] :
  myTime.append(str(i))
newFolder = 'cfgs_{TIME}'.format(TIME='_'.join(myTime))
os.system('mkdir {TIME}'.format(TIME=newFolder))
os.system('cp setup_cff.py {TIME}'.format(TIME=newFolder))
os.system('mkdir {TIME}/{OUT}'.format(TIME=newFolder, OUT=options.outdir))
os.chdir(os.getcwd()+'/' + newFolder)
try :
  os.path.isfile(options.input)
except :
  print 'unknown input filelist'
  sys.exit(0)
flist   = open(options.input)
infiles = flist.readlines()
sample  = '_'.join([infiles[0].split("/")[2] , infiles[0].split("/")[3] , infiles[0].split("/")[5]]) 

myFirstFile = int(options.start)
if myFirstFile > 1:
  infiles = infiles[int(myFirstFile)-1:]

if int(njobs) == -1 : 
  njobs = len(infiles)/int(nfiles) + 1*(len(infiles)%int(nfiles)>0) 
  
for j in range(njobs): 
  k = myFirstFile + j
  print " k is " + str(k)
  f   = open(options.cfg)
  f1  = open('{M}_{NUM}.py'.format(M=options.cfg.replace(".py","").replace("../","").rstrip(), NUM=str(k)), "w")
  for line in f:
      newline = None
      newfile = None
      if 'inputFileNameSub' in line:
          line = line.replace('inputFileNameSub', '').rstrip()
          for i in infiles[k*nfiles:(k+1)*nfiles]:
            if i=="": continue
            print >> f1,'\t\t"'+ i.rstrip() + '",'
      if 'outputFileNameSub' in line:
          newline = line.replace('outputFileNameSub', '"{TIME}/{OUT}/{TRIG}_{TAG}_{NUM}.root"'.format(TIME=newFolder, OUT=options.outdir, NUM=str(k), TAG=sample, TRIG=options.hlt )).rstrip()
#           newline = line.replace('outputFileNameSub', '"root://eoscms//eos/cms/store/user/fiorendi/{OUTDIR}/{TRIG}_{TAG}_{NUM}.root"'.format(OUTDIR=options.outdir, NUM=str(k), TAG=sample, TRIG=options.hlt )).rstrip()
      if newline:
        print >> f1,newline.strip()
      if not (newline or newfile):
        print >> f1,line.rstrip() 
  
  sh   = open("../script.sh")
  shName = "script_{NUM}.sh".format(NUM=str(k))
  sh1  = open(shName,"w")
  os.system("chmod +x {SH}".format(SH=shName)) 
  for shline in sh:
      mycfg = None
      myres = None
      if 'hlt.py'  in shline:
          shline = shline.replace('cfgs/hlt.py', newFolder + '/{M}_{NUM}.py'.format(M=options.cfg.replace("../","").replace(".py","").rstrip(), NUM=str(k))).rstrip()
      if 'hlt.log' in shline:
          shline = shline.replace('hlt.log', newFolder + '/{M}_{NUM}.log'.format(M="log", NUM=str(k))).rstrip()
      if 'RESFILE' in shline:
          myres  = shline.replace('RESFILE', '{TIME}/{OUT}/{TRIG}_{TAG}_{NUM}.root'.format(TIME=newFolder, OUT=options.outdir, NUM=str(k), TAG=sample, TRIG=options.hlt )).replace('EOSDIR', '/eos/cms/store/user/fiorendi/{OUTDIR}/'.format(OUTDIR=options.outdir)).rstrip()
#           myres  = shline.replace('RESFILE', '{TIME}/{OUT}/{TRIG}_{TAG}_{NUM}.root'.format(TIME=newFolder, OUT=options.outdir, NUM=str(k), TAG=sample, TRIG=options.hlt )).replace('EOSDIR', 'root://eoscms//eos/cms/store/user/fiorendi/{OUTDIR}/{TRIG}_{TAG}_{NUM}.root'.format(OUTDIR=options.outdir, NUM=str(k), TAG=sample, TRIG=options.hlt )).rstrip()
      if myres:
        print >> sh1, myres.rstrip()
#       if shline:
#         print >> sh1,shline.strip()
      else: 
        print >> sh1, shline.rstrip()
      
  
  os.system("bsub -q {QUEUE} {SH}".format(QUEUE=options.queue, SH=shName))    
#   os.system("chmod 755 root://eoscms//eos/cms/store/user/fiorendi/{OUTDIR}/{TRIG}_{TAG}_{NUM}.root".format(OUTDIR=options.outdir, NUM=str(k), TAG=sample, TRIG=options.hlt))