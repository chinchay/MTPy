import os
import time
import sys
import subprocess
from subprocess import Popen, PIPE
import os.path
from os import path
# You need `pot_blank_ternary.mtp`, `cycles.py`,  `to_relax.cfg`, `jobtrain.sh`, `train.cfg` (create it by `touch train.cfg` if it's your first cycle)
# You can run this code with:
# salloc --time=10:00:00 --ntasks=1 --cpus-per-task=1 --mem-per-cpu=2G --account=def-rmelnik
# module --force purge && module load StdEnv/2016.4 && module load nixpkgs/16.09 intel/2019.3 intelmpi/2019.3.199 && module load python/3.6.3 && source /home/chinchay/projects/def-rmelnik/chinchay/mydocs/venvs/jupyter_py3/bin/activate
# python -u cycles.py | tee mylog.txt ##(to save screen messages into mylog.txt)

####################################################################
def sleepWhileTrainJobIsPD():
    count = 0
    while not path.exists("2_myTraining/training.txt"):
        print("2_myTraining/training.txt does not exist. mlp train job is still in queue or it has not been sent!")
        time.sleep(60)
        count += 1
        # if it has been sleeping for 1 hour:
        if count > 60:
            sys.exit()
        #
    #
#

def getTrainJobId():
    command = 'squeue -u chinchay | grep "jobtrain.sh"'
    jobId = ''
    count = 0
    while (jobId == '') and (count < 20):
        try:
            line = os.popen(command).read()
            if line == '':
                jobId = 'nojobfound'
            else:
                jobId = line.split()[0]
            #
        except:
            print("something went wrong when asking slurm for jobId..., I will try again")
            sleep(60)
            count += 1
        #
    #
    if count >= 20:
        print("error while asking slurm for jobId, stopping...")
        sys.exit()
    #
    return jobId
#

def sleepWhileDFTJobsRunning():
    command = 'squeue -u chinchay | grep "META" | wc -l'
    nSqueue = -1
    count = 0
    while (nSqueue == -1) and (count < 20):
        try:
            nSqueue = int(os.popen(command).read().split()[0])
        except:
            print("something went wrong when asking slurm for DFT jobs..., I will try again")
            sleep(30)
            count += 1
        #
    #
    if count >= 20:
        print("error while asking slurm for DFT jobs, stopping...")
        sys.exit()
    #
    if nSqueue >= 1:
        return True
    #
    return False
#

def trainingFinished():
    endTrainingLine = "_______________________________________________"
    command = "tail -2 2_myTraining/training.txt"
    checkLine = os.popen(command).read().split()[0]
    return ( endTrainingLine == checkLine )
#

def getLastLineTraining():
    return os.popen("tail -1 2_myTraining/training.txt").read().split("\n")[0]
#

def isTrainingStuck(checkTrainTime=15, process=None):    
    prevLine = getLastLineTraining()
    time.sleep(checkTrainTime)
    nextLine = getLastLineTraining()
    print("prevLine...........: " + prevLine)
    print("nextLine...........: " + nextLine)
    print("")
    return (prevLine == nextLine)
#

def checkErrorFile():
    g = open("2_myTraining/errorsByPythonPopen.txt", "r")
    checkLine = g.readline().split(" \n")[0]
    g.close()
    return checkLine
#

def hasBFGSascendingError():
    g = open("2_myTraining/errorsByPythonPopen.txt", "r")
    checkLine = g.readline().split(" \n")[0]
    g.close()
    #
    errorBFGSaccending = "ERROR: BFGS: stepping in accend direction detected."
    return ( checkErrorFile() == errorBFGSaccending )
#

def killTrainJob():
    jobID = getTrainJobId()
    if jobID != 'nojobfound': # there is no training job!
        command = "scancel " + jobID
        os.system(command)
    #
#

def train(checkTrainTime):
    # clean in case you are coming here from a previous stuck train: there will be a training.txt!
    command = "cd 2_myTraining/  && " +\
              "rm -f errorsByPythonPopen.txt temp1.cfg training.txt state.mvs selected.cfg"
    os.system(command)

    command = "wc -l 2_myTraining/train.cfg"
    nLines = int(os.popen(command).read().split()[0])

    if nLines > 3000:
        command = "cd 2_myTraining/  && " +\
                  "sbatch jobtrain.sh"
        os.system(command)

        sleepWhileTrainJobIsPD() # <<== here it waits while training.txt is not found!
        # Now TrainJob is running

        time.sleep(checkTrainTime) # <<== let training start
        while not trainingFinished():
            isStuck = isTrainingStuck(checkTrainTime)
            if isStuck:
                killTrainJob()
                return isStuck, None
            #
        #
        # Training has finished, and isStuck=False
        return isStuck, False  # returns: isStuck, hasBFGSascendingError()
        # ********** I SHOULD IMPROVE THE PREVIOUS LINE CODE *****
        # ********** IT SHOULD BE return isStuck,  hasBFGSascendingError(), but I have not implemente  hasBFGSascendingError() for slurm job sendin :/ only for process in the next part 
        #
    else:
        f = open("2_myTraining/training.txt", "w")
        g = open("2_myTraining/errorsByPythonPopen.txt", "w")
        p = Popen(["srun", "mlp", "train", "2_myTraining/pot.mtp", "2_myTraining/train.cfg"], stdout=f, stderr=g)

        time.sleep(checkTrainTime) # <<== let training start
        while not trainingFinished():
            isStuck = isTrainingStuck(checkTrainTime, p)
            if isStuck:
                p.kill()
                print("I had to kill the training process because it is stuck")
                f.close()
                g.close()                
                return isStuck, None
            #
        #
        # Training has finished, and isStuck=False
        f.close()
        g.close()
        # Training has finished, and isStuck=False
        return isStuck, hasBFGSascendingError()  # returns: isStuck, hasBFGSascendingError()
    #
#
####################################################################

def initialize():
    os.system("mkdir 2_myTraining/")
    os.system("mkdir 4_toRelax/")
    os.system("mkdir 5_afterActiveLearning/")

    # command = "cd 5_afterActiveLearning/  && " +\
    #           "cp -r /home/chinchay/projects/def-rmelnik/chinchay/mydocs/META/ .  && " +\
    #           "cp -r /home/chinchay/projects/def-rmelnik/chinchay/mydocs/saved/ . && " +\
    #           "cd META/  &&  yes | clean.run  &&  cp ../saved/* ."  ## << it doesn't work :/
    command = "cd 5_afterActiveLearning/  && " +\
              "cp -r /home/chinchay/projects/def-rmelnik/chinchay/mydocs/META/ .  && " +\
              "cp -r /home/chinchay/projects/def-rmelnik/chinchay/mydocs/saved/ . && " +\
              "cd META/  &&  cp ../saved/* ."
    os.system(command)

    command = "cd 4_toRelax/  && " +\
              "cp /home/chinchay/projects/def-rmelnik/chinchay/mydocs/paper1/e_trainedforcarbonwire/4_toRelax/relax.ini .  && " +\
              "cp ../to_relax.cfg ."
    os.system(command)

    # command = "cd 2_myTraining/  && " +\
    #           "cp /home/chinchay/projects/def-rmelnik/chinchay/mydocs/paper1/agnr/1/2_myTraining/pot_blank_binary.mtp .  && " +\
    #           "cp pot_blank_binary.mtp pot.mtp  && " +\
    #           "cp ../train.cfg .  && " +\
    #           "mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg"
    # os.system(command)

    # command = "cd 4_toRelax/  && " +\
    #           "cp ../2_myTraining/pot.mtp .  && " +\
    #           "cp ../2_myTraining/state.mvs .  && " +\
    #           "mlp relax relax.ini --cfg-filename=to_relax.cfg --min-dist=0.5  && " +\
    #           "cat selected.cfg_*  > selected.cfg"
    # os.system(command)

    existsPrevious = False
    if path.exists("previousPot.mtp"):
        print("previousPot.mtp found")
        os.system("cp previousPot.mtp 2_myTraining/pot.mtp")

        if path.exists("previousState.mvs") and path.exists("previousTrain.cfg"):
            print("previousState.mvs and previousTrain.cfg found")
            os.system("cp previousState.mvs 2_myTraining/state.mvs")
            os.system("cp previousTrain.cfg 2_myTraining/train.cfg")
            existsPrevious = True
        #
        else:
            print("previousState.mvs or previousTrain.cfg not found, introduce it in workingDirectory")
            sys.exit()
        #
    else:
        os.system("cp pot_blank_ternary.mtp 2_myTraining/pot_blank_ternary.mtp")
        os.system("cp 2_myTraining/pot_blank_ternary.mtp 2_myTraining/pot.mtp")
        os.system("cp train.cfg 2_myTraining/train.cfg")
    #

    os.system("cp jobtrain.sh 2_myTraining/jobtrain.sh")
    #
    return existsPrevious
#

def selectionStep():
    # clean from previous selection step:
    command = "cd 2_myTraining/  && " +\
              "rm -f diff.cfg POSCAR*"
    os.system(command)

    command = "cd 2_myTraining/  && " +\
              "cp ../4_toRelax/selected.cfg .  && " +\
              "srun mlp select-add pot.mtp train.cfg selected.cfg diff.cfg  && " +\
              "srun mlp convert-cfg diff.cfg POSCAR --output-format=vasp-poscar"
    os.system(command)
    print("Finished Selection step")
    checkTostop()
#

def sendJobs(nJobs):
    command = "cd 5_afterActiveLearning/META/  && " +\
              "submit.run " + str(nJobs)
    os.system(command)
    
    continueSleeping = True
    while continueSleeping:
        time.sleep(120)
        # continueSleeping = shouldContinueSleeping()
        continueSleeping = sleepWhileDFTJobsRunning()
    #
    ###
    command = "cd 5_afterActiveLearning/META/  && " +\
              "resubmit.run " + str(nJobs)
    string   = "No failed/unfinished jobs; nothing to resubmit"
    st2check = os.popen(command).read().split("\n")[0]
    if st2check != string:
        continueSleeping = True
        while continueSleeping:
            time.sleep(120)
            # continueSleeping = shouldContinueSleeping()
            continueSleeping = sleepWhileDFTJobsRunning()
        #    
    ######

    # The following did not work because I cannot get the correct signal for waiting until finished jobs
    # I tried to implement the advice here https://blog.miguelgrinberg.com/post/how-to-make-python-wait
    # other advices here https://hpc-discourse.usc.edu/t/how-to-wait-until-the-job-is-finished/242/2
    # command = 'cd 5_afterActiveLearning/META/  && ' +\
    #           'echo "submit.run ' + str(2) + ' --wait &" > myexecute.sh  && ' +\
    #           'echo "wait" >> myexecute.sh  && ' +\
    #           'bash myexecute.sh'
    # os.system(command)

    # ########
    # # wait

    # command = 'cd 5_afterActiveLearning/META/  && ' +\
    #           'echo "resubmit.run ' + str(2) + ' --wait &" > myexecute.sh  && ' +\
    #           'echo "wait" >> myexecute.sh  && ' +\
    #           'bash myexecute.sh'
    # os.system(command)

    ########
    # wait
    # if "nothing to resubmit" continue
#

def getNumberJobs2send(nPoscars):
    # this is a temporary solution...
    return (nPoscars // 6) + 1 # job_script.sh have 1 hour of wall time, and each sub-job takes ~10min to finish (for 9atoms)
#

def dftStep():
    command = "cd 5_afterActiveLearning  && rm -rf META/"
    os.system(command)

    command = "cd 5_afterActiveLearning/  && " +\
              "cp -r /home/chinchay/projects/def-rmelnik/chinchay/mydocs/META/ .  && " +\
              "cd META/  &&  cp ../saved/* .  && " +\
              "cp ../../2_myTraining/POSCAR* . "
    os.system(command)

    command = "ls 5_afterActiveLearning/META/POSCAR* | wc -l"
    nPoscars = int(os.popen(command).read().split()[0])

    os.system("rm -f 5_afterActiveLearning/META/table.dat")
    f = open("5_afterActiveLearning/META/table.dat", "a")
    for i in range(nPoscars):
        command = "cd 5_afterActiveLearning/META/  &&  python vasp2qe.py -i \"POSCAR" + str(i) + "\""
        os.system(command)
        
        line = "srun pw.x < ../scf_" + str(i + 1) + ".in > scf_" + str(i + 1) + ".out\n"
        f.write(line)
    #
    f.close()
    #
    print("sending jobs now...")
    assert nPoscars > 0
    nJobs = getNumberJobs2send(nPoscars)
    sendJobs(nJobs)
    print("It seems I have finished jobs")
    print("Finished DFT step")
    checkTostop()
#

def copyFromDFT2Training():
    # clean from previous training step:
    command = "cd 2_myTraining/  && " +\
              "rm -f errorsByPythonPopen.txt temp1.cfg training.txt state.mvs selected.cfg"
    os.system(command)

    command = "cd 5_afterActiveLearning/META/  && " +\
              "python scriptQE2cfg.py  && " +\
              "cd ../../2_myTraining/  && " +\
              "rm -f POSCAR* && rm -f diff.cfg && rm -f temp1.cfg && rm -f selected.cfg && rm -f training.txt  && " +\
              "cp ../5_afterActiveLearning/META/train.cfg train2.cfg  && " +\
              "cat train2.cfg >> train.cfg  && " +\
              "rm train2.cfg"
    os.system(command)

    command = "cat 5_afterActiveLearning/META/scf_*  >> qeins.txt"
    os.system(command)

    command = "cat 5_afterActiveLearning/META/RUN*/scf_*  >> qeouts.txt"
    os.system(command)

    print("files from dft to train ready")
    checkTostop()
#

def trainingStep(checkTrainTime):
    print("mlp train ...")
    #
    # errorBFGSaccending = "ERROR: BFGS: stepping in accend direction detected."
    isTrainingOK = False
    for i in range(5):        
        if not isTrainingOK:
            gotStuck, BFGSaccendingError = train(checkTrainTime)
            if gotStuck:
                # mlp got stuck, runs infinitely. `training.txt` has stuck in the same line
                isTrainingOK = False
                print("Training got stuck")
            #
            elif BFGSaccendingError:
                isTrainingOK = False
                print("There is a BGFS ascending error :/")
            else:
                isTrainingOK = True
                print("it seems there are no errors?...")
            #
        #
    #
    if not isTrainingOK:
        print("trying another trick...")
        # let's try another trick:
        # if the same error appear again, take a fresh pot.mtp:
        command = "cp 2_myTraining/pot_blank_ternary.mtp 2_myTraining/pot.mtp"
        os.system(command)
        gotStuck, BFGSaccendingError = train(checkTrainTime)
        if not gotStuck:
            if BFGSaccendingError:
                #? You should increase configs in 4_toRelax/to_reala.cfg, and start all again
                print("error errorBFGSaccending again?")
                print("You should increase configs in 4_toRelax/to_relax.cfg, and start all again")
                print("I am going to stop here... :/")
                sys.exit()
        else:
            print("Training keeps getting stuck")
            print("You should increase configs in 4_toRelax/to_relax.cfg, and start all again")
            print("I am going to stop here... :/")
            sys.exit()
        #
    #
    #
    print("Finished Training step")
    checkTostop()
#

def updateTrainedPotential():
    if path.exists("Trained.mtp_"):
        os.system("ls >> mylog.txt")
        os.system("mv Trained.mtp_  2_myTraining/pot.mtp")# <<-- mtp was run in an outside directory, so this corrects the location file
        os.system("echo ''")
        os.system("ls >> mylog.txt")
    elif path.exists("2_myTraining/Trained.mtp_"):
        os.system("ls 2_myTraining/  >> mylog.txt")
        os.system("mv 2_myTraining/Trained.mtp_  2_myTraining/pot.mtp")
        os.system("echo ''")
        os.system("ls 2_myTraining/  >> mylog.txt")
    else:
        print("I cannot find Trained.mtp_ , stopping...")
        sys.exit()
    #

    command = "cat 2_myTraining/training.txt >> 2_myTraining/alltrainings.txt"
    os.system(command)

    print("Finished updating trained potential")
    checkTostop()
#

def calcGrade():
    command = "cd 2_myTraining/  && " +\
              "srun mlp calc-grade pot.mtp train.cfg train.cfg temp1.cfg"
    os.system(command)
    #
    print("Finished calc-grade")
    checkTostop()
#

def relaxStep():
    command = "cd 4_toRelax/  && " +\
              "rm -f select*  && " +\
              "cp ../2_myTraining/pot.mtp .  && " +\
              "cp ../2_myTraining/state.mvs .  && " +\
              "srun mlp relax relax.ini --cfg-filename=to_relax.cfg --min-dist=0.5  && " +\
              "cat selected.cfg_*  > selected.cfg"
    os.system(command)

    print("Finished Relaxation step")
    checkTostop()
#

def checkTostop():
    if path.exists("stop.txt"):
        print("I found 'stop.txt', so I should stop...")
        sys.exit()
    #
#

# os.system("mlj")##??

def countCfgsFile(fileName):
    command = 'grep "BEGIN_CFG" ' + fileName + ' | wc -l'
    return int(os.popen(command).read().split()[0])
#

def shouldContinue(fileName):
    return ( countCfgsFile(fileName) != 0 )
#

## Test later the interface below. Needs improvement, and initialization!
# from itertools import cycle
# def mtprun(selection, maxNcycles, checkTrainTime):
#     steps = ["selectionStep()", "dftStep()", "copyFromDFT2Training()", "trainingStep(checkTrainTime)", "updateTrainedPotential()", "calcGrade()", "relaxStep()"]
#     pool  = cycle(steps)

#     l = range(len(steps))
#     pool  = cycle(l)
#     for _ in range(selection - 1):
#         s = next(pool)
#     #

#     continuar = True
#     for i in range(maxNcycles):
#         funStr = steps[next(pool)]
#         if continuar:
#             while funStr != "relaxStep()":
#                 print( "step: " + funStr )
#                 eval(funStr)
#                 funStr = steps[next(pool)]
#             #
#             eval(funStr) # = "relaxStep()"
#             continuar = shouldContinue("4_toRelax/selected.cfg")
#             print("end of loop i = ", i )
#         #
#     #
#     print("finish loop")
# #
# selection = 1
# maxNcycles = 2
# checkTrainTime = 15
# mtprun(selection, maxNcycles, checkTrainTime)





maxNcycles = 15
#nJobs = 1
checkTrainTime = 30

if not path.exists("notinit.txt"):
    existsPrevious = initialize()
    if existsPrevious:
        print("I am going to relax because I found previous trained potential ")
        relaxStep()
        print("finished initialization using previous potential")
    else:
        nCfgsTrain = countCfgsFile("2_myTraining/train.cfg")
        print("nCfgsTrain = ", nCfgsTrain)
        if nCfgsTrain != 0:
            trainingStep(checkTrainTime)
            updateTrainedPotential()
        #
        calcGrade()
        relaxStep()
        print("finished initialization2")
    #
#
continuar = True
for i in range(maxNcycles):
    if continuar:
        selectionStep()
        dftStep()
        copyFromDFT2Training()
        trainingStep(checkTrainTime)
        updateTrainedPotential()
        calcGrade()
        relaxStep()
        continuar = shouldContinue("4_toRelax/selected.cfg")
        print("end of loop i = " + str(i) )
#     #
print("finish loop")
# #
