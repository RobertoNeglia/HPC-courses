# Bash Essentials

Docker gives you the possibility to build and interact with another OS inside your machine. Running the command `docker exec -it hpc-courses /bin/bash` will connect you with the bash shell of an OS where you will find all the tools for the hpc courses.

## Shell
Shell is a command line interface for running programs on a computer. The user types a bunch of commands at the prompt, the shell runs the programs for the user and then displays the output. The commands can be either directly entered by the user or read from a file called the shell script or shell program.

## Exploring the Terminal
Whenever we open a terminal window, we see our last login credentials and a Shell prompt. The Shell prompt appears whenever the shell is ready to accept the input. It appears as `jellyfish@561c4c540a77 ~ #`, where `jellyfish` is your username, `561c4c540a77` is the name of the machine, and `~` is your current directory.

## Getting started
To get a little hang of the bash, let’s try a few simple commands:
- `echo`: returns whatever you type at the shell prompt similar to `print` in Python, or `printf` / `disp` in Matlab.
- `date`: displays the current time and date.
- `clear`: clean the terminal.

## Variables and Environmental Variables
You can assign a value to a variable with the equal sign **(no spaces!)**, for instance type `A=1`. You can then retrieve its value using the dollar sign and curly braces, for instance to display it the user may type `echo ${A}`. Some variables can affect the way running processes will behave on a computer, these are called environmental variables. For this reason, some variables are set by default, for instance to display the user home directory type `echo ${HOME}`. To set an environmental variable just prepend `export`, for instance `export PATH="/usr/sbin:$PATH"` adds the folder `/usr/sbin` to the `PATH` environment variable. `PATH` specifies a set of directories where executable programs are located.

## Basic Bash Commands
- `pwd` stands for **Print working directory** and it points to the current working directory, that is, the directory that the shell is currently looking at. It’s also the default place where the shell commands will look for data files.
- `ls` stands for a **List** and it lists the contents of a directory. ls usually starts out looking at our home directory. This means if we print ls by itself, it will always print the contents of the current directory.
- `cd` stands for **Change directory** and changes the active directory to the path specified.
- `mkdir` stands for **Make directory** and is used to make a new directory or a folder.
- `mv` stands for **Move** and it moves one or more files or directories from one place to another. We need to specify what we want to move, i.e., the source and where we want to move them, i.e., the destination.
- `touch` command is used to create new, empty files. It is also used to change the timestamps on existing files and directories.
- `rm` stands for **Remove** and it removes files or directories. By default, it does not remove directories, but if used as `rm -r *` within a directory, then every directory and file inside that directory is deleted (`*` is a special characters that matches everything).

### Exercises
**Suggestion:** You can check the documentation for each command with the flag `--help`. For instance, `ls --help` will list you all the options for listing files in a directory. Please, try to rely only on the documentation to solve the exercises.

- Go to your home folder (*Suggestion:* you can either use `~` or `$HOME`)
- Create a folder named `test1`
- Go inside `test1` and create a directory `test2`
- Go inside `test2` and then up one directory (*Suggestion:* `..` indicates the upper directory)
- Create the following files `file1.txt`, `file2.txt`, `file3.dat`, `file4.md`, `README.md`, `.hidden.txt`
- List all files in the directory
- List also the hidden files in the directory
- List only files with txt extension (*Suggestion:* use `*` wildcard)
- List only files that start with `file`
- List files with `1`, `2`, `3` or `4` in the name (*Suggestion:* use `[1-4]` wildcard)
- Move the `README.md` in `test2`
- Move all txt files in `test2` in one command
- Remove `file3.dat`
- Remove all contents of `test2` and the folder itself in one commands

## Download and Uncompress a Matrix
With `wget` you can retrieve content from web servers. For instance, you can download a matrix from the matrix market with `wget https://math.nist.gov/pub/MatrixMarket2/NEP/mhd/mhd416a.mtx.gz`.
To unzip the file, simply type `gzip -dk mhd416a.mtx.gz`

## Other Commands, Pipelines and Filters
- `cat` stands for **Concatenate** and it reads a file and outputs its content. It can read any number of files, and hence the name concatenate.
- `wc` is short for **Word count**. It reads a list of files and generates one or more of the following statistics: newline count, word count, and byte count.
- `grep` stands for **Global regular expression print**. It searches for lines with a given string or looks for a pattern in a given input stream.

We can add operators between commands in order to chain them togheter.
- The pipe operator `|` (vertical bar), is a way to send the output of one command as an input to another command. E.g. with `cat /etc/passwd | grep jellyfish` you can check system informations about your user.
- The redirect operator `>` is a way to send the output of one command to a file (which might not exist). E.g. with `ls > files-in-this-folder.txt` you save a file with the list of files in this directory.
- The append operator `>>` is a way to append the output of one command to a file (which might not exist).

### Exercises
- Create a file with the current date (one command) and display its content
- Count the number of lines in the matrix `mhd416a.mtx` (*Suggestion:* use `cat`, `wc` and `|`)
- List the entries of the matrix that are smaller than 1e-10 in absolute value. You can assume that all values are in exponential format and all values are greater than 1e-100 in absolute value. Count how many entries satisfy this criteria. (*Suggestion:* use `cat`, `grep`, `wc` and `|` )

