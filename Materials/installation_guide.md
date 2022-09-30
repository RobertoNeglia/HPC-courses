# Installation of dependencies: mk modules

## 1. What are mk modules?

[mk modules](https://github.com/elauksap/mk) bundle a set of scientific libraries compiled under the same toolchain. Once installed, they provide the command module, that has several subcommands:

```
module load <module name> 
```

loads the requested module. This creates a set of environment variables storing relevant paths for that library (e.g. `mkEigenPrefix`, `mkEigenInc`, ...). Use

- `export | grep mk`  to obtain a list
- `module list`: to show a list of currently loaded modules
- `module avail`: to show a list of all available modules (loaded or not)
- `module --help`: to show a list of all the commands

## 2. Installation
![Installation Flowchart](./assets/installation-flowchart.png)

### 2.1. Install Docker
As a first step, install the [Docker container environment](https://www.docker.com/). You can follow the instruction on the [official guide](https://docs.docker.com/get-docker/). Please read it thoroughly.

**IMPORTANT (Windows users):** Docker is avaiable on Windows only for the following versions:
* Windows 11 64-bit: Home or Pro version 21H2 or higher, or Enterprise or Education version 21H2 or higher.
* Windows 10 64-bit: Home or Pro 21H1 (build 19043) or higher, or Enterprise or Education 20H2 (build 19042) or higher.

To check you version press WINDOWS + R, enter `winver` and press OK. If your version is not compatible with Docker, you will need to use **dual boot**. See Section 2.4.

### 2.2. Pull the Docker image
From a terminal with admin privileges, run the command

```
docker pull elauksap/hpc_courses
```

The image is just a snapshot of the state of a Linux OS, it is like a saving point from where you want to start. You can check your images with `docker image ls`.

### 2.3. Use the Docker image 
To use your image you need to create a Docker container. To make a parallel with virtual machines, the Docker image is like the .iso of the OS, but then you have to install it. We want to create a container with the image we have just downloaded, give it a name (`--name hpc-courses`) to remember its function and share a folder with the host so that we can exchange file easily (`-v /path/to/host/folder:/home/jellyfish/shared-folder`). The complete command is:

```
docker run --name hpc-courses -v /path/to/host/folder:/home/jellyfish/shared-folder -it elauksap/hpc_courses
```

**WARNING:** to avoid problems `/path/to/host/folder` should not contain white spaces or special characters. For instance you can make your shared folder with the command `mkdir shared-folder` and than `/path/to/host/folder` would be `C:/Users/matteo/shared-folder` on Windows or `~/shared-folder` on Linux-like OS.

You have now created a container. To turn on the container type:

```
docker start hpc-courses
```
To enter into the container run:

```
docker exec -it hpc-courses /bin/bash
```
You can leave the container and return to your OS with `exit`. You can check your containers and their status with the command

```
docker ps -a
```
If the status of the container is `Up`, you can stop it with

```
docker stop hpc-courses
```
Once you have created your container remember to **do not** use again the commad `run` but just `start`. Otherwise you will create every time a new container. If you want to remove a container you creaded for mistake you can run:

```
docker rm <name-of-the-container>
```

Always remember that documentation is your best friend! Do not panic, just type:
```
docker --help
```

### 2.4 Dual Boot
Dual-booting is the act of installing two operating systems on a single computer, and being able to choose which one to boot. You will need to install a Linux distribution, we suggest Ubuntu 20.04 LTS or 22.04 LTS. 

Before proceeding we suggest you to backup you data. To install a Dual Boot you can follow the [official guide](https://help.ubuntu.com/community/WindowsDualBoot) (which is a bit dated) or this [unofficial tutorial](https://itsfoss.com/install-ubuntu-1404-dual-boot-mode-windows-8-81-uefi/).

After you're done, you can install Docker on the Linux OS, following steps 2.1 to 2.3 above.

## 3. Editing files

1. Install a text editor on the host OS (i.e. the one in which you installed Docker). We advise you to use [Visual Studio Code](https://code.visualstudio.com/), which is available for Linux, Windows and MacOS.

2. Create and edit files inside the shared folder, so that they will be visible to both the host and the Docker container.

3. You will edit files from the host OS, and compile them from inside the Docker container.

## 4. Test the installation

1. Using VS Code, open the shared folder and create a file `test-installation.cpp` with content:

```
#include <Eigen/Eigen>
#include <iostream>

int main(int argc, char** argv)
{
        std::cout << "Successfully included Eigen." << std::endl;
        return 0;
}
```

2. Change the current directory to the shared folder `cd /home/jellyfish/shared-folder`.

3. In the container, make sure the Eigen module is loaded: `module load eigen`.

3. Compile and run the test:

```
g++ -I ${mkEigenInc} test-installation.cpp -o test-installation
./test-installation
```
