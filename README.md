# Setting up CUT-RUN Pipeline on Riviera
The code in the CUT-RUN repository is intended to be run on HPC clusters, such as Riviera. This document discusses how to set up the repository, Conda environments, etc. after getting access to Riviera.

## Table of contents
- [Using Riviera](#using-riviera)
	- [Logging in with MobaXterm](#logging-in-to-riviera-with-mobaxterm)
	- [Logging in from terminal/the command line](#logging-in-to-riviera-from-the-command-line)
	- [Transferring files from your computer (without MobaXterm)](#transferring-files-without-mobaxterm)
	- [Transferring files to/from the NAS](#transferring-files-directly-between-the-nas-and-riviera)
- [Downloading the CUT&RUN GitHub Files](#downloading-the-cut-run-github-files)
- [Setting up Conda](#setting-up-conda-miniconda3)
	- [Using conda](#using-conda)
	- [Creating the CUT&RUN conda environments](#creating-the-cut-run-conda-environments)
	- [Removing conda environments](#removing-conda-environments)
	- [Creating your own conda environment](#creating-your-own-conda-environment)
	- [Adding/removing packages in conda](#installing-packages-to-your-conda-environment)
- [Code run on login: .bashrc](#code-run-on-login-bashrc)
	- [Auto-starting conda](#auto-starting-conda)
	- [Setting environment variables](#setting-environment-variables)
	- [Loading modules](#loading-modules)
	- [Custom commands](#creating-custom-commands)
- [How to use Slurm](#how-to-use-slurm)

## Using Riviera
**YOU MUST USE THE GLOBALPROTECT VPN WHEN CONNECTING TO RIVIERA IF YOU ARE NOT ON CAMPUS.**

If you have not already done so, you must contact DSRI to get an account. Fill out the form on the page for Riviera here: https://www.research.colostate.edu/dsri/hpc-riviera/
If it has been several days and you haven't already gotten a notice about your account, you may need to contact Bill Carpenter (contact information is on the website) to get your login credentials.

Once you have your credentials, you can log on to Riviera. If you are using Windows, MobaXterm (https://mobaxterm.mobatek.net/) is recommended for its GUI file explorer and ease of downloading/uploading files. There aren't many great alternatives for Mac or Linux, but you can run it using the Wine compatibility tool. (go to https://www.winehq.org/ if you are on Mac, or install `wine` with your package manager on Linux. Then, download and run the portable MobaXterm .exe)

### Logging in to Riviera with MobaXterm
Click the button labeled "Session" in the top menu. A new window will pop up; select the SSH option at the top. Enter `riviera.colostate.edu` in the "Remote host" box, check "Specify username", and enter your given username in that box. Click the "OK" button, and you will see a new shell session appear.
Double click on the shell session; you will see a login prompt on the terminal. Enter your password for Riviera. (it will not show up as you type it)

### Logging in to Riviera from the command line
To log into Riviera without a tool like MobaXterm, run the following command from the command line, replacing `username` with your username:
``` bash
ssh username@riviera.colostate.edu
```

You will be prompted for your password. It will not show up as you type, so type carefully.

### Transferring files without MobaXterm
If you cannot use MobaXterm or another GUI-based SSH tool, you can upload and download files with the `scp` command. Here is a template for uploading files to Riviera from a MacOS or Linux machine. Because it uses the `-r` flag, it can be used to copy individual files or directories. Wildcards (`*`) are allowed, but can occasionally cause issues. The tilde `~` is used to denote the home directory, but you can use an absolute path by starting with a `/` if needed. Replace the paths and USERNAME with the correct values.
``` bash
scp -r path/to/file/or/directory USERNAME@riviera.colostate.edu:~/path/to/destination
```

Here is a similar template for downloading files from Riviera:
``` bash
scp -r USERNAME@riviera.colostate.edu:~/path/to/file/or/directory path/to/destination
```

If you have directories or files with spaces in their names (which should generally be avoided), then try wrapping the paths in quotes like so:
``` bash
scp -r "path/to/file or directory" USERNAME@riviera.colostate.edu:"~/path/to your/destination"
```

### Transferring files "directly" between the NAS and Riviera
Sometimes it is inconvenient to transfer files through your own computer. (this is what happens when you drag files directly from the NAS to a window of MobaXterm)  It can be useful to set up an asynchronous process that handles these file transfers. This tends to be faster, and it allows you to transfer files without keeping your computer on. This requires access to the AIDLNGS01 server, or any other server that has a *mounted* NAS drive. The following command, which should be run on Riviera, uses the AIDLNGS01 server to transfer files from the NAS to Riviera:

``` bash
sshpass -p "YourAIDLNGS01Password" scp -r AidlngsUsername@AIDLNGS01.cvmbs.colostate.edu:/home/lab_data/rosenberg_lab/NAS_DIRECTORY/path/to/files path/to/destination
```

There are Slurm versions of this command saved on the NAS in `NAS_Server_Shared/shared_scripts/ZM_custom`. I would suggest putting them in your home directory on Riviera, so that you can call them with `sbatch ~/transfer_from_async.sh` or `sbatch ~/transfer_to_async.sh` (the `~` refers to your home directory)

## Downloading the CUT-RUN GitHub files
We will use `git` to download the CUT&RUN pipeline files straight from the GitHub.
1. `cd` into the directory where you want to install the CUT&RUN pipeline.
```bash
# 2. Load the git module
module load git
# 3. Clone the CUT-RUN github into your current directory
git clone https://github.com/CRosenbergCode/CUT-RUN.git
```

## Setting up Conda (miniconda3)
Conda is a package and environment manager, allowing you to switch between different "environments" that have different programs installed. This is useful when your programs may conflict with one another.

We are following the instructions listed on https://www.anaconda.com/docs/getting-started/miniconda/install#macos-linux-installation.

``` bash
# 1. Make sure you're in your home directory
cd ~
# 2. Download miniconda3 installation script
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
# 3. Run the installation script
bash ~/Miniconda3-latest-Linux-x86_64.sh
```
4. Scroll through Anaconda's TOS with the Return key, then enter `yes` when prompted to agree to the TOS.
5. Press Return to use the default installation directory (`~/miniconda3`)
6. Enter `yes` to enable auto-initialization. (you can change this later in the [.bashrc section](#code-run-on-login-bashrc)).
7. The first time you do this, conda will not auto-initialize. You can make it initialize by re-running your .bashrc script:
``` bash
source ~/.bashrc
```
8. Finally, you have to configure conda, mainly for installing packages correctly.
``` bash
# Adding a common conda channel, conda-forge.
conda config --add channels conda-forge
# If a package exists in other channels, we want to install it from there instead of from Bioconda
conda config --append channels bioconda
# Recommended package
conda-forge install tree
```

### Using conda
The main conda commands you will use are `conda activate` and `conda deactivate`.

1. When you need to switch into a conda environment, first check that you're in the base environment by looking to the left of your username in the terminal. It should say `(base)`. If it says something else, you must run the command `conda deactivate`. If there is nothing in the parentheses to the left of your name, then conda is not running, and you must start it with the command `conda activate`.

2. Once conda is running, you can change into a new environment by using `conda activate EnvironmentName`, replacing "EnvironmentName" with the proper name of your environment. If you are unsure what its name is, run `conda env list` to see a list of available environments.
```bash
conda activate EnvironmentName
```

### Creating the CUT-RUN conda environments
All the necessary conda environments have already been created as .yml files in the CUT-RUN files. Go to `CUT-RUN-main/CondaEnvs/` and run the following command for each .yml file you need to create an environment from. You will then be prompted to provide permission to install different packages, respond with `y`.
``` bash
conda env create --file YourCondaEnv.yml
```

### Removing conda environments
If you have made a mistake and need to remove a conda environment, run `conda env list` to ensure you have the right name, then delete the environment by name with the following command.. **MAKE SURE YOU HAVE THE ENVIRONMENT DEACTIVATED BEFORE REMOVING IT.**
``` bash
# 1. Check that you have the right environment name
conda env list
```
``` bash
# 2. Remove the environment
conda remove -n YourEnvironmentName --all
```

### Creating your own conda environment
There are many options for creating conda environments. If you are at the point at which you're creating your own conda environments, you are probably comfortable reading the [conda docs page on managing environments](https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html). Otherwise, the most basic way to create an environment is using the following command.
``` bash
conda create --name MyEnvironmentName
```

You can then activate this environment with the following conda command:
``` bash
conda activate MyEnvironmentName
```

### Installing packages to your conda environment
Conda packages can be installed with the command `conda install`. To find out the exact name of a package, you should always search it up. For example, let's walk through installing the `rmats` package to our current environment:
1. First, find the rmats package online. Search "conda rmats". This brings up a bioconda page that shows up the correct command to install rmats: `conda install bioconda::rmats` (Notice that it is part of the `bioconda` channel. You will often install packages from here)
2. On Riviera, ensure you are in the environment you want to install packages in. Check the [Using conda](#using-conda) instructions above. In this case, let's say we are in the "MyEnv" environment:
``` bash
conda activate MyEnv
```
3. Install the rmats package:
``` bash
conda install bioconda::rmats
```
4. You will be prompted to install the necessary packages. Enter `y` to install.
5. The packages will install in the MyEnv environment!

As a rule of thumb, do not install too many unrelated programs in the same environment. They may conflict and you will have issues installing new packages, and you may find yourself in a dependency situation that [looks like this](https://xkcd.com/1987/).

### Removing packages from your conda environment
If you need a list of available packages in the current environment, run `conda list`. Once you have your package name, remove it with the following command. You will be prompted to remove the packages. Enter `y` to remove them:
``` bash
conda remove MyPackage
```
You can remove multiple packages by separating the package names with spaces:
``` bash
conda remove MyPackage MyOtherPackage
```

## Code run on login: .bashrc
Hidden files and directories are preceded by a dot. The `.bashrc` file is a hidden file in your home directory that contains code to be run on startup. This can be a lot of fun to edit and personalize, (for example, you can have something that displays welcome messages or tells you the weather for the day) but we'll be focusing on setting environment variables and how to enable/disable auto-starting conda.

To edit your .bashrc file, run the following command:
``` bash
nano ~/.bashrc
```
You will enter a command line-based text editor called nano. When you are finished making edits, type Ctrl+X and then enter `Y` to save your edits, or `N` to not.

<img src="https://i.imgur.com/Y3cDw16.png" width="500"/>

*This user has a .bashrc script that draws a picture on login*


### Auto-starting conda
If you previously did not opt to auto-start conda on log-in, you will need to add the code for it to your .bashrc. Copy and paste the following code into your .bashrc, replacing `zmikol` with your own Riviera username. If you would like to turn off auto-initialization, delete this code from your .bashrc.
``` bash
# >>> conda initialize >>>  
# !! Contents within this block are managed by 'conda init' !!  
__conda_setup="$('/nfs/home/zmikol/miniconda3/bin/conda'  'shell.bash'  'hook'  2> /dev/null)"  
if [  $?  -eq  0  ];  then  
	eval  "$__conda_setup"  
else  
	if [  -f  "/nfs/home/zmikol/miniconda3/etc/profile.d/conda.sh"  ];  then  
		. "/nfs/home/zmikol/miniconda3/etc/profile.d/conda.sh"  
	else  
		export PATH="/nfs/home/zmikol/miniconda3/bin:$PATH"  
	fi  
fi  
unset __conda_setup  
# <<< conda initialize <<<
```

### Setting environment variables
To set an environment variable on startup, use the following code:
- `export myvar=myvalue`, where myvar is your variable, and myvalue is the value you are setting it to.

The most helpful environment variable to set on Riviera is the TMPDIR variable. Many programs will default to using Riviera's `/tmp` directory, which has a disappointing 2 gigabytes of allocated space. You should create your own temp directory. You can name it anything you like, but the following is suggested: (this creates a "hidden" directory in your home folder)
``` bash
mkdir ~/.temp
```

Then, you can **add the following code to your .bashrc**, replacing `username` with your own Riviera username.
``` bash
export TMPDIR="/nfs/home/username/.temp"
```

If your temp directory has a different name or is elsewhere, `cd` into it, then copy the output of the `pwd` command, and set that as the value of your TMPDIR variable in .bashrc.

### Loading modules
Certain tools on Riviera must be loaded as modules, with the most important one being Slurm. To load a module, you can use the `module load module_name` command. (where module_name is the name of your module) Add the following code to your .bashrc:
``` bash
module load slurm
module load git
```

You can check out the other available modules by running `module spider`, and you can search for specific modules, such as CUDA, like so: `module spider cuda`. Use the up/down arrow keys to navigate, and `q` to quit.

### Creating custom commands
This section is optional, but can help create shortcuts that save time. You can create a custom command in bash using the following pattern:
``` bash
function myfunction() {
	# code goes here
}
```

By putting custom commands in your .bashrc, you can call them from the command line simply by entering their name. In this case, we would be able to run the `myfunction` command, although it wouldn't do anything, since we haven't added any code inside it.

Some custom commands may just be shortcuts to calling other scripts. For example, the [asynchronous file transfer scripts](#transferring-files-directly-between-the-nas-and-riviera) above could be turned into easy commands like so. `$1` and `$2` refer to the first and second arguments sent to the command.
``` bash
function getnas() {
	sbatch ~/transfer_from_async.sh $1 $2
}

function sendnas() {
	sbatch ~/transfer_to_async.sh $1 $2
}
```

Another useful command is a quick slurm job scheduler for running big jobs that don't have their own Slurm script. (For example, copying or unzipping a 300GB directory)   This requires you to add the `salias_helper.sh` script (located in `NAS_Server_Shared/shared_scripts/ZM_custom`) to your home directory first, and rename it to `.salias_helper.sh` (adding a period at the start so it's out of sight)
``` bash
function salias() {
	sbatch ~/.salias_helper.sh $@
}
```
Then, you can run your big commands, like unzipping a massive file, as slurm jobs like so: `salias unzip my_massive_file.zip`

## How to use Slurm
Slurm is a jobs manager for HPC systems like Riviera. There is fairly comprehensive documentation you can read at https://slurm.schedmd.com/documentation.html. This will be a brief overview of the most common commands, `sbatch`, `srun`, `scancel`, `sacct`, `squeue`, and `sinfo`, as well as how to write a script for Slurm.
- **sbatch -** This is the most commonly-used command to run a slurm job. You must provide it with the file to be run, as well as any job parameters like dependencies. The basic pattern is `sbatch [options] script [args]`. 
	- Most options can be written in the header of a script, which will be shown later. Some options must be determined when the job is run, like a dependency. An example of this can be seen in the array example scripts in the CUT-RUN github.
	- After the options, you write the path to the script you will be executing as a Slurm job.
	- After the script, you can provide the arguments that will be passed it. For example, the `03_SICC.sh` script requires a path to the bam files to be converted to bigwigs. It must be run like so: `sbatch 03_sicc.sh path/to/bams/folder`
- **srun -** srun is an interactive version of sbatch, allowing the user to see the program as it executes and provide inputs as needed. This command is less common in general, and entirely unused in the CUT-RUN pipeline, but it is still important to know. It uses the same syntax as sbatch.
- **scancel -** scancel stops a slurm job based on ID. You may need this if you realize you have made an error. You can find a job ID by running `sacct` or `squeue`. Job IDs are also printed to the terminal upon being scheduled. Although it has other options as shown on the docs, you will likely only ever need to know the following two uses:
	- `scancel 12345` to cancel job 12345.
	- `scancel {31415..31425}` to cancel jobs 31415 to 31425 (Only if you have scheduled several jobs at once and they all have sequential job IDs. Double-check to make sure you aren't trying to cancel someone else's job!)
- **sacct -** sacct shows your recent slurm jobs and their status. (running, failed, completed, waiting, etc.) In most cases, simply entering `sacct` to the terminal on its own is enough. However, there are some cases where you may want to see more/less than its regular output.
	- If you've been scheduling a lot of jobs recently, you may just want to subset it to the recent ones. In that case, you can run it with `tail`, using twice the number of lines as recent jobs you want to see. `sacct | tail -n10`
	- If you're checking on a job that was run a day prior, you'll need to specify that with the `--starttime` option, and possibly the `--endtime` option as well. Formatting of days is `YYYY-MM-DD`, and in order to use specific times of the day you must use the pattern `YYYY-MM-DDTHH:MM:SS`. Seconds may be omitted. That can be a bit confusing, so here is an example of a user checking the status of jobs that ran after *March 14th, 2025 at 3:32PM* but before *April 29th, 2025 at 2:50AM and 29 seconds*:
		- `sacct --starttime 2025-03-14T15:32 --endtime 2025-04-29T02:50:29` (note the use of 24-hour time)
- **squeue -** squeue shows current slurm jobs from all users, including those currently waiting on dependencies. It can simply be run by running `squeue` in the terminal. You can also combine it with `grep` to find jobs specifically belonging to you, like so: `squeue | grep yourusername`
- **sinfo -** Generally not used, but sinfo can be helpful when scheduling jobs to see which partitions are currently being taken up most. If, for example, the week-long-highmem partition is being taken up and you have a program that normally uses that partition, you can check on the status of other partitions by running `sinfo`. If week-long-cpu is open and your program can still run with reduced available memory, you can run it on that partition instead.
