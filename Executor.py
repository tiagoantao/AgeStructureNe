#Implements asynchronous executors for local machines and condor
#
#The local executor is only partially async, if there is too much load
#   it will block!

import multiprocessing
import os
import subprocess
import time


class Condor:
    ''' Class for Condor.
    '''
    def __init__(self):
        ''' Constructor

            universe is vanilla OR standard ...
                Only vanilla is implemented (but standard should work)
        '''
        self.universe = "vanilla"
        self.job = 0

    def java(self, command, parameters, requirements):
        '''  Submits a condor java job.
        '''
        self.job += 1
        template = '''
universe=java
requirements={0!s}
should_transfer_files=yes
when_to_transfer_output=on_exit
transfer_input_files={1!s}
log=condor{2:d}.log
error=condor{3:d}.error
output={4!s}
executable={5!s}
arguments={6!s}
java_vm_args={7!s}
queue
        '''.format(requirements, ",".join(self.send), self.job, self.job, self.out, command, parameters, self.javaVMArgs)

        rootName = "condor{0:d}".format(self.job)
        jobName = rootName + ".job"
        logName = rootName + ".log"
        w = open(jobName, "w")
        w.write(template)
        w.close()
        subprocess.call("condor_submit {0!s}".format(jobName), shell=True)
        self.logName = logName

    def submit(self, command, parameters):
        '''  Submits a condor job.
        '''
        self.job += 1
        template = '''
universe={0!s}
requirements=(Arch=="X86_64" || Arch=="INTEL") && ( OpSys=="LINUX")
rank = kflops
should_transfer_files=yes
when_to_transfer_output=on_exit
transfer_input_files={1!s}
log=condor{2:d}.log
error=condor{3:d}.error
output={4!s}
executable={5!s}
arguments={6!s}
queue
        '''.format(self.universe, ",".join(self.send), self.job, self.job, self.out, command, parameters)

        rootName = "condor{0:d}".format(self.job)
        jobName = rootName + ".job"
        logName = rootName + ".log"
        w = open(jobName, "w")
        w.write(template)
        w.close()
        subprocess.call("condor_submit {0!s}".format(jobName), shell=True)
        self.logName = logName


class Local:
    ''' The local executor.

        This executor will block if there are no more slots available!
    '''
    def __init__(self, limit):
        '''

             if limit is an int>0 then it is the expected load average NOT
                 to be used, for instance if there are 32 cores and there
                 is a limit of 6, the system will try to never ago above 26.
                 A float between 0 and 1 will be interpreted as a the
                 fraction of CPUs to be used, e.g., with 32 cores, 0.25
                 will use at most 8
                 A negative value will be interpreted as the maximum number
                 of processes that can be executed in parallel.

        '''
        self.limit = limit
        self.cpus = multiprocessing.cpu_count()
        self.running = []
        pass

    def cleanDone(self):
        '''Removes dead processes from the running list.
        '''
        myDels = []
        for rIdx, p in enumerate(self.running):
            if p.poll() is not None:
                myDels.append(rIdx)
        myDels.reverse()
        for myDel in myDels:
            del self.running[myDel]

    def wait(self, forAll=False):
        '''Blocks if there are no slots available

           forAll - Also waits if there is *ANY* job running (i.e.
                    block/barrier)
        '''
        self.cleanDone()
        numWaits = 0
        if self.limit > 0 and type(self.limit) == int:
            cond = "len(self.running) >= self.cpus - self.limit"
        elif self.limit < 0:
            cond = "len(self.running) >= - self.limit"
        else:
            cond = "len(self.running) >= self.cpus * self.limit"
        while eval(cond) or (forAll and len(self.running) > 0):
            time.sleep(1)
            self.cleanDone()
            numWaits += 1

    def submit(self, command, parameters):
        '''Submits a job
        '''
        self.wait()
        if hasattr(self, "out"):
            out = self.out
        else:
            out = "/dev/null"
        if hasattr(self, "err"):
            err = self.err
        else:
            err = "/dev/null"
        if err == "stderr":
            errSt = ""
        else:
            errSt = "2> " + err
        p = subprocess.Popen("{0!s} {1!s} > {2!s} {3!s}".format(command, parameters, out, errSt),
                             shell=True)
        self.running.append(p)
        if hasattr(self, "out"):
            del self.out
        if hasattr(self, "err"):
            del self.err


class Pseudo:
    ''' The pseudo executor.

        This executor will Dump of list of nohup nice commands
    '''
    def __init__(self, outFile="/tmp/pseudo{0:d}".format((os.getpid()))):
        '''
           outFile is where the text is written
        '''
        self.outFile = open(outFile, "a")
        pass

    def submit(self, command, parameters):
        '''Submits a job
        '''
        self.outFile.write("nohup /usr/bin/nice -n19 {0!s} {1!s} > {2!s}\n".format(command, parameters, self.out))
        self.outFile.flush()

    def __del__(self):
        self.outFile.close()


class LSF:
    ''' The LSF executor.

    '''
    def __init__(self):
        '''Constructor

        '''
        self.running = []
        self.queue = "normal"  # Default queue name is "normal"
        self.mem = 4000  # Request 4GB as a default
        self.numPasses = 0
        self.outDir = os.path.expanduser("~/tmp")
        self.cnt = 1

    def cleanDone(self):
        '''Removes dead processes from the running list.
        '''
        ongoing = []
        statusFile = "/tmp/farm-{0:d}".format((os.getpid()))
        os.system("bjobs > {0!s} 2>/dev/null".format(statusFile))
        f = open(statusFile)
        f.readline()  # header
        for l in f:
            toks = filter(lambda x: x != "", l.rstrip().split(" "))
            ongoing.append(int(toks[0]))
        os.remove(statusFile)
        myDels = []
        for rIdx, p in enumerate(self.running):
            if p not in ongoing:
                myDels.append(rIdx)
        myDels.reverse()
        for myDel in myDels:
            del self.running[myDel]

    def wait(self, forAll=False, beCareful=60):
        '''Blocks according to some condition

           forAll - Also waits if there is *ANY* job running (i.e.
                    block/barrier)

           beCareful - Wait X secs before starting. This is because
                       tasks take time to go into the pool.
        '''
        time.sleep(beCareful)
        if forAll:
            while len(self.running) > 0:
                self.cleanDone()
                time.sleep(1)

    def submit(self, command, parameters="", myDir=os.getcwd()):
        '''Submits a job
        '''
        M = self.mem * 1000
        job = "bsub -G malaria-dk -P malaria-dk -q %s "
        job += "-o quickrun.%s.out -e quickrun.%s.err "
        job += "-J quickrun.%s -M %d -R "
        job += "'select[type==X86_64 && mem>%d] "
        job += "rusage[mem=%d]' \"cd %s ; %s %s \""
        job = job % (self.queue, self.cnt, self.cnt, self.cnt, M,
                     self.mem, self.mem, myDir, command, parameters)

        statusFile = "/tmp/farm-{0:d}.{1:d}".format(os.getpid(), self.cnt)
        os.system(job + " >{0!s} 2>/dev/null".format(statusFile))
        f = open(statusFile)
        l = f.readline()
        job = int(l[l.find("<") + 1:l.find(">")])
        f.close()
        os.remove(statusFile)
        self.cnt += 1

        self.running.append(job)
        self.numPasses += 1


class SGE:
    ''' The SGE executor.

    '''
    def __init__(self, mailUser=None):
        '''Constructor

        '''
        self.running = []
        self.queue = "normal"  # Default queue name is "normal"
        self.mem = 1000  # Request 1GB as a default
        self.outDir = os.path.expanduser("~/tmp")
        self.cnt = 1
        self.project = "anopheles"
        self.outDir = "/tmp"
        self.mailOptions = "a"
        self.mailUser = mailUser
        self.maxProc = 1000
        self.hosts = []

    def cleanDone(self):
        '''Removes dead processes from the running list.
        '''
        ongoing = []
        statusFile = "/tmp/farm-{0:d}".format((os.getpid()))
        os.system("qstat > {0!s} 2>/dev/null".format(statusFile))
        f = open(statusFile)
        f.readline()  # header
        f.readline()  # header
        for l in f:
            toks = filter(lambda x: x != "", l.rstrip().split(" "))
            ongoing.append(int(toks[0]))
        os.remove(statusFile)
        myDels = []
        for rIdx, p in enumerate(self.running):
            if p not in ongoing:
                myDels.append(rIdx)
        myDels.reverse()
        for myDel in myDels:
            del self.running[myDel]

    def wait(self, forAll=False, beCareful=60):
        '''Blocks according to some condition

           forAll - Also waits if there is *ANY* job running (i.e.
                    block/barrier)

           beCareful - Wait X secs before starting. This is because
                       tasks take time to go into the pool.
        '''
        time.sleep(beCareful)
        self.cleanDone()
        if forAll:
            while len(self.running) > 0:
                self.cleanDone()
                time.sleep(1)

    def submit(self, command, parameters="", myDir=os.getcwd()):
        '''Submits a job
        '''
        if self.mailUser is not None:
            mail = "-m {0!s} -M {1!s}".format(self.mailOptions, self.mailUser)
        while len(self.running) > self.maxProc:
            self.wait(beCareful=5)
        hosts = ""
        if len(self.hosts) > 0:
            hosts = " -q "
        for host in self.hosts:
            hosts += "\*@{0!s}".format(host)
            if host != self.hosts[-1]:
                hosts += ","
        job = "qsub {0!s} {1!s} -S /bin/bash -V -P {2!s} -cwd -l h_vmem={3:d}m {4!s} {5!s}".format(
            mail, hosts, self.project, self.mem, command, parameters)
        statusFile = "/tmp/farm-{0:d}.{1:d}".format(os.getpid(), self.cnt)
        os.system(job + " >{0!s} 2>/dev/null".format(statusFile))
        f = open(statusFile)
        l = f.readline()
        job = int(l.split(" ")[2])
        f.close()
        os.remove(statusFile)
        self.cnt += 1

        self.running.append(job)
