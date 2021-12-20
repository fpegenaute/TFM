import os
import pwd
import subprocess
import time
from bin.config import AF2config, SLURMconfig, SSHconfig

def get_filename_ext(filepath):
    """
    Given a file path, split it by dots and get the extension as the last element
    Also, take the basename in the path and split it by the point, get the first
    element as the name of the file
    """

    extension = filepath.split(".")[-1]
    filename = os.path.basename(filepath).split('.')[0]
    
    return filename, extension

def number_of_jobs_in_queue(squeue="squeue"):
    """
    This functions returns the number of jobs in queue for a given
    user.
    """

    # Initialize #
    user_name = "$USER"

    process = subprocess.check_output([squeue, "-u", user_name])

    return len([line for line in process.split("\n") if user_name in line])


def submit_command_to_queue(command, queue="normal", max_jobs_in_queue=None, queue_file=None, dummy_dir=".", submit="sbatch", squeue="squeue"):
    """
    This function submits any {command} to a cluster {queue}.
    @input:
    command {string}
    queue {string} by default it submits to any queue (partition)
    max_jobs_in_queue {int} limits the number of jobs in queue
    queue_file is a file with information specific of the cluster for running a queue
    """
    import hashlib

    if max_jobs_in_queue is not None:
        while number_of_jobs_in_queue(squeue) >= max_jobs_in_queue: time.sleep(5)

    cwd = os.path.join(dummy_dir)
    if not os.path.exists(cwd): 
        os.makedirs(cwd)

    script= os.path.join(cwd,"submit_"+hashlib.sha224(command).hexdigest()+".sh")

    # Use the config for SLURM given in a file and add your command
    if queue_file is not None:
      fd=open(script,"w")

      with open(queue_file,"r") as queue_standard:
        data=queue_standard.read()
        fd.write(data)
        fd.write("%s\n\n"%(command))
      fd.close()
      queue_standard.close()
      
      # choose which workload manager to use
      if queue is not None:
       if  submit=="qsub":
          os.system("%s -q %s %s" % (submit, queue,script))
       elif submit=="sbatch":
          os.system("%s -p %s %s" % (submit, queue,script))
       else:
          os.system("%s %s"% (submit,script))
      else:
        os.system("%s %s"% (submit,script))
    else:
      if queue is not None:
        os.system("echo \"%s\" | %s -q %s" % (command, submit, queue))
      else:
        os.system("echo \"%s\" | %s" % (submit,command))

def submit_AF_to_SLURM(query_fasta, outdir, workload_manager="sbatch", dummy_dir=".", max_jobs_in_queue=None ):
    """
    This function submits AF2 command to a cluster via {workload_manager}.
    
    @input:
    query_fasta: fasta file to be analyzed
    outdir: output for *AlphaFold*, not SLURM
    workload_manager: SLURM by default
    dummy_dir=default directory for SLURM to look at, defaul: the current
    max_jobs_in_queue {int} limits the number of jobs in queue

    """
    import hashlib
    from datetime import date

    today = str(date.today())
    



    if max_jobs_in_queue is not None:
        while number_of_jobs_in_queue("squeue") >= max_jobs_in_queue: time.sleep(5)

    cwd = os.path.join(dummy_dir)
    if not os.path.exists(cwd): 
        os.makedirs(cwd)

    script= os.path.join(cwd,"alpha_"+hashlib.sha224("Alphafold").hexdigest()+".sh")

    command = f"""bash $alphafold_path/run_alphafold.sh -d {AF2config["AF2datadir"]} 
                                                        -o {outdir} 
                                                        -m {AF2config["AF2preset_monomer"]} 
                                                        -f {query_fasta} 
                                                        -t {today} 
                                                       --is_prokaryote_list={AF2config["AF2_prokaryote"]} 
                                                        -g {AF2config["AF2_useGPU"]}
    """ 

    # Use the config for SLURM given in a file and add your command
    with open(script,"w") as batch_script:
        batch_script.write(SLURMconfig)
        batch_script.write("%s\n\n"%(command))
        batch_script.write("conda deactivate\n")

    # Connect with the HPC Cluster:
    # command = "ssh -T marvin"
    # result = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stdin=subprocess.PIPE).communicate()
    # Better like this: https://kb.iu.edu/d/aews, https://stackoverflow.com/questions/8382847/how-to-ssh-connect-through-python-paramiko-with-ppk-public-key
      
    os.system("%s %s" % (workload_manager,script))

