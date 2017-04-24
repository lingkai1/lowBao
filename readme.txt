1.Overall Introduction
(1) Use the script for compilation and packing. If the compilation fails, solve this problem by yourself.
(2) If the compilation succeeds, the executable binary file cdn will be generated in the bin/ directory.
(3) Invoke and debug your program using the ./cdn /xxx/topo.txt /xxx/result.txt command, where topo.txt is the input file and result.txt is the output file.
(4) Submit the cdn.tar.gz package on the contest website after the successful debugging, and query your score later.
the answer structure like this:
cdn.tar.gz/ (end with tar.gz or zip, the file name can not contain spaces and characters.)
|-- bin/ 
|-- build/ 
|-- cdn/  
|-- build.sh	   Can not be modified or deleted
|-- readme.txt

2.SDK Directory Structure
SDK-gcc.zip/
|-- bin/                         Directory storing the executable binary file. The shell script deletes this directory before compilation
|                                and then creates it again. Therefore, script execution is not affected by this directory.
|-- build/                       Build directory. The shell script deletes this directory before compilation and then creates it again.
|                                Therefore, script execution is not affected by this directory.
|-- cdn/                         Code directory.
|     |-- lib/                   lib header file directory, which cannot be modified. No files can be added to this directory.    
|     |     |-- lib_io.h         Header file for reading and writing files, which cannot be modified.   
|     |     |-- lib_time.h       Header file for printing time, which cannot be modified.   
|     |-- CMakeLists.txt         cmake, which cannot be modified.   
|     |-- cdn.cpp                main function source file, which cannot be modified.  
|     |-- io.cpp                 Source file providing file reading and writing functions, which cannot be modified. 
|     |-- deploy.cpp             Source file whose code is to be written.
|     |-- deploy.h               Header file whose code is to be written.
|-- build.sh                     Script for compilation and packing, which cannot be modified.
|-- readme.txt                   The file you are reading.

3.Description of build.sh Shell Script
 You can execute this script for compilation and packing. If the compilation is correct, the binary file cdn will be generated in the bin/ directory, and a compressed code source file (cdn.tar.gz) will be generated in the the same directory. 
NOTE:
 (1) The shell script will delete the bin/ and build/ directories as well as all files and sub-directories under these two directories. Do not save any files to these two directories. 
 (2) If you want to use the shell script function, ensure the integrity of contents in the SDK-gcc/ directory. Do not change the name of any directory and file, and do not move any directory or file.

4.(Optional) Manual Operation Description (This part can be ignored if you use the shell script) 
 (1) Create the build_private/ directory under the SDK-gcc/ directory, and edit the makefile file in the build_private/ directory. 
 (2) Go to the build_private/ directory and run the make command.
 (3) Pack and compress the code source file to cdn.tar.gz.
NOTE: 
 (1) Do not store your makefile file in the build/ directory. This file will be deleted once batch.sh is invoked. 
 (2) Ensure that all source files can be seen once the package is decompressed. Otherwise, the compilation fails. 

5.SDK Code Description
 Now, necessary information has been explained. What you need to do is as follows:  
 (1) Implement the XXXX interface in the deploy.cpp file.
 (2) You need to write the result by writing_result method in the format required by the contest title.
 (3) If there is no solution, the output is NA.
 The SDK provides the functions of reading files, writing files in required format, and printing the start and end time. 
 NOTE: The function of reading files refers to reading the link information files by lines to the memory. The file storage format remains to be the character string format. The storage format is related to the algorithm design, which will spark your thoughts. 

6.Note:
 Please note that only the code source file you have modified or added can be submitted in the preliminary. The compilation will be performed on the contest question server. Therefore, special attentions must be paid to the following things: 
 (1) The code development must be based on the SDK. Otherwise, the compilation fails. 
 (2) In the SDK source files, only deploy.cpp and deploy.h files can be modified. Do not modify other files, otherwise it may compile failed;
 (3) New files must be added to the same directory of the deploy.cpp file.
