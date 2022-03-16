#compile and run stdin1 from command line providing file ending in .cpp
#stdin 2 is the optional input fastq if compiling and running with input text
#stdin 3 is the optional output fastq if compiling and running with input text
cpp_file=$1
argv_1=$2
argv_2=$3

echo 'compiling' $cpp_file
#/usr/bin/g++ $cpp_file -o ${cpp_file%%.cpp} && ./${cpp_file%%.cpp} $argv_1
/usr/bin/g++ -std=c++11 $cpp_file -o ${cpp_file%%.cpp} && ./${cpp_file%%.cpp} $argv_1 $argv_2
